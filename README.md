R1CS Distributifier
=======

Simple tool to prepare R1CS for more efficient distributed proof generation

## Instructions for Hekaton
1. Run
```bash
r1cs-distributify prepare circuit.r1cs circuit.hgr
```
to prepare a hypergraph.
2. Use KaHyPar (or maybe hMetis) to partition the hypergraph. An example configuration would be
```bash
KaHyPar -o cut -m direct -p ../../../config/cut_kKaHyPar_sea20.ini -e 0.05 -k 4 -h circuit.hgr -w true
```
3. Run
```bash
r1cs-distributify finalize -c circuit.r1cs -w circuit.json -n 4 -p circuit.hgr.part4.epsilon0.05.seed-1.KaHyPar -o circuit
```
to partition the single R1CS into multiple files. 3 files will be produced per block (`r1cs`, `json` witness, and `meta`).
