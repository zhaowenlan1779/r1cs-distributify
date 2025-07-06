use ark_bn254::Fr;
use ark_ff::Field;
use circom_compat::{read_binary_wtns, read_witness, write_witness, R1CSFile};
use clap::{Args, Parser, Subcommand};
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Args)]
struct PrepareArgs {
    /// R1CS circuit file (e.g. circuit.r1cs)
    input: String,

    /// Output hypergraph file (e.g. circuit.hgr)
    output: String,
}

#[derive(Args)]
struct FinalizeArgs {
    /// R1CS circuit file (e.g. circuit.r1cs)
    #[arg(short, long)]
    circuit: String,

    /// R1CS witness file (e.g. circuit.json or circuit.wtns)
    #[arg(short, long)]
    witness: String,

    #[arg(short, long)]
    n_blocks: usize,

    /// Partition file from KaHyPar (e.g. circuit.hgr.part4.epsilon0.05.seed-1.KaHyPar)
    #[arg(short, long)]
    partition: String,

    /// Path of the common prefix of the output files (e.g. circuit)
    #[arg(short, long)]
    output: String,
}

#[derive(Subcommand)]
enum Commands {
    Prepare(PrepareArgs),
    Finalize(FinalizeArgs),
}

use std::io::Write;
use std::iter::zip;

fn prepare(cli: PrepareArgs) {
    let reader = BufReader::new(File::open(cli.input).unwrap());
    let r1cs = R1CSFile::<Fr>::new(reader).unwrap();

    let mut file = File::create(cli.output).unwrap();
    let mut variable_appearances = vec![vec![]; r1cs.header.n_wires as usize];
    for (i, (a, b, c)) in r1cs.constraints.iter().enumerate() {
        let mut indices = a
            .iter()
            .chain(b.iter())
            .chain(c.iter())
            .map(|(var, _)| *var)
            .collect::<Vec<_>>();
        indices.sort();
        indices.dedup();
        for index in indices {
            variable_appearances[index].push(i);
        }
    }
    let mut count = 0;
    for var in variable_appearances.iter().skip(1) {
        if var.len() > 1 {
            count += 1;
        }
    }
    writeln!(&mut file, "{} {}", count, r1cs.header.n_constraints).unwrap();
    for var in variable_appearances.iter().skip(1) {
        if var.len() > 1 {
            for i in var {
                write!(&mut file, "{} ", i + 1).unwrap();
            }
            writeln!(&mut file, "").unwrap();
        }
    }
}

fn check_r1cs_satisfied(r1cs: &R1CSFile<Fr>) {
    assert!(r1cs.constraints.par_iter().all(|(a, b, c)| {
        let a_val = a
            .iter()
            .map(|(var, coeff)| *coeff * r1cs.witness[*var])
            .sum::<Fr>();
        let b_val = b
            .iter()
            .map(|(var, coeff)| *coeff * r1cs.witness[*var])
            .sum::<Fr>();
        let c_val = c
            .iter()
            .map(|(var, coeff)| *coeff * r1cs.witness[*var])
            .sum::<Fr>();
        a_val * b_val == c_val
    }));
}

fn finalize(cli: FinalizeArgs) {
    let reader = BufReader::new(File::open(cli.circuit).unwrap());
    let mut r1cs = R1CSFile::<Fr>::new(reader).unwrap();

    let is_json_witness = cli.witness.to_ascii_uppercase().ends_with(".JSON");
    let witness_reader = BufReader::new(File::open(cli.witness).unwrap());

    r1cs.witness = if is_json_witness {
        read_witness::<Fr>(witness_reader)
    } else {
        read_binary_wtns::<Fr>(witness_reader).unwrap()
    };

    let partition_reader = BufReader::new(File::open(cli.partition).unwrap());
    let mut partition = partition_reader
        .lines()
        .map(|line| {
            let line = line.unwrap();
            line.parse::<usize>().unwrap()
        })
        .collect::<Vec<_>>();
    assert_eq!(partition.len(), r1cs.constraints.len());

    // First reorder the blocks, so the one with the fewest constraints come first
    let mut variable_uses = vec![BTreeSet::new(); r1cs.header.n_wires as usize];
    let n_blocks = cli.n_blocks;

    let mut constraint_count = vec![0; n_blocks];
    for block in &partition {
        constraint_count[*block] += 1;
    }
    let mut block_constraint_count = constraint_count.into_iter().enumerate().collect::<Vec<_>>();
    block_constraint_count.sort_by_key(|(a, b)| *b);
    let mut new_block_map = vec![0; n_blocks];
    for (i, (block, _)) in block_constraint_count.iter().enumerate() {
        new_block_map[*block] = i;
    }

    for val in &mut partition {
        *val = new_block_map[*val];
    }

    for ((a, b, c), block) in zip(&r1cs.constraints, &partition) {
        for (var, _) in a {
            variable_uses[*var].insert(*block);
        }
        for (var, _) in b {
            variable_uses[*var].insert(*block);
        }
        for (var, _) in c {
            variable_uses[*var].insert(*block);
        }
    }

    let mut variable_assignments = vec![vec![]; n_blocks];
    for (var, uses) in variable_uses.iter().enumerate().skip(1) {
        if uses.len() == 1 {
            let block = uses.iter().next().unwrap();
            variable_assignments[*block].push(var);
        }
    }
    let unique_variables_count = variable_assignments
        .iter()
        .map(|x| x.len())
        .collect::<Vec<_>>();
    let mut borrowed_variables = vec![vec![]; n_blocks];
    for (var, uses) in variable_uses.iter().enumerate().skip(1) {
        if uses.len() == 1 {
            continue;
        }
        let block = if uses.len() == 0 {
            (0..n_blocks)
                .min_by_key(|block| {
                    variable_assignments[*block].len() + borrowed_variables[*block].len()
                })
                .unwrap()
        } else {
            // Need to assign to first subcircuit
            *uses.iter().next().unwrap()
        };
        variable_assignments[block].push(var);
        if uses.len() > 1 {
            for other_block in uses {
                if *other_block != block {
                    borrowed_variables[*other_block].push(var);
                }
            }
        }
    }
    let mut reorder_map = vec![HashMap::new(); n_blocks];
    for block in 0..n_blocks {
        reorder_map[block].insert(0, 0);
        for (new_var, old_var) in variable_assignments[block].iter().enumerate() {
            reorder_map[block].insert(*old_var, new_var + 1);
        }
        let num_owned_vars = variable_assignments[block].len() + 1;
        for (borrow_idx, old_var) in borrowed_variables[block].iter().enumerate() {
            reorder_map[block].insert(*old_var, num_owned_vars + borrow_idx);
        }
    }

    let mut sub_r1cs = (0..n_blocks)
        .map(|_| R1CSFile::<Fr>::new_bn254())
        .collect::<Vec<_>>();
    for block in 0..n_blocks {
        sub_r1cs[block].header.n_wires =
            (variable_assignments[block].len() + borrowed_variables[block].len() + 1) as u32;
        sub_r1cs[block].witness.push(Fr::ONE);
        for var in variable_assignments[block]
            .iter()
            .chain(borrowed_variables[block].iter())
        {
            sub_r1cs[block].witness.push(r1cs.witness[*var].clone());
        }
    }
    for ((a, b, c), block) in zip(&r1cs.constraints, &partition) {
        let a = a
            .iter()
            .map(|(var, coeff)| (*reorder_map[*block].get(var).unwrap(), coeff.clone()))
            .collect::<Vec<_>>();
        let b = b
            .iter()
            .map(|(var, coeff)| (*reorder_map[*block].get(var).unwrap(), coeff.clone()))
            .collect::<Vec<_>>();
        let c = c
            .iter()
            .map(|(var, coeff)| (*reorder_map[*block].get(var).unwrap(), coeff.clone()))
            .collect::<Vec<_>>();
        sub_r1cs[*block].constraints.push((a, b, c));
    }
    for block in 0..n_blocks {
        sub_r1cs[block].header.n_constraints = sub_r1cs[block].constraints.len() as u32;
    }

    for block in 0..n_blocks {
        check_r1cs_satisfied(&sub_r1cs[block]);

        let r1cs_out =
            BufWriter::new(File::create(format!("{}.{}.r1cs", cli.output, block)).unwrap());
        sub_r1cs[block].write(r1cs_out).unwrap();

        let witness_out =
            BufWriter::new(File::create(format!("{}.{}.json", cli.output, block)).unwrap());
        write_witness(&sub_r1cs[block].witness, witness_out).unwrap();

        let mut meta_out = File::create(format!("{}.{}.meta", cli.output, block)).unwrap();
        writeln!(
            &mut meta_out,
            "{} {} {}",
            unique_variables_count[block] + 1,
            variable_assignments[block].len() - unique_variables_count[block],
            borrowed_variables[block].len()
        )
        .unwrap();
        for idx in unique_variables_count[block]..variable_assignments[block].len() {
            writeln!(&mut meta_out, "{}", variable_assignments[block][idx]).unwrap();
        }
        for var in &borrowed_variables[block] {
            writeln!(&mut meta_out, "{}", *var).unwrap();
        }
    }
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::Prepare(args) => prepare(args),
        Commands::Finalize(args) => finalize(args),
    }
}
