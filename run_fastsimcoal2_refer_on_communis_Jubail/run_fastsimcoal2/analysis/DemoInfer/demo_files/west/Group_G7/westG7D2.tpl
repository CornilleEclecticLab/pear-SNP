//Parameters  for the coalescence simulation program : simcoal.exe
4 samples to simulate :
//Population effective sizes (number of genes)
Ncauc
NcomD
NcomP
Npyra
//Samples sizes, samples ages, and inbreeding levels
42
28
22
18
//Growth rates  : negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
4
//Migration matrix 0
0 0 M02 M03
0 0 0 0
M20 0 0 0
M30 0 0 0
//Migration matrix 1
0 0 0 M03
0 0 0 0
0 0 0 0
M03 0 0 0
//Migration matrix 2
0 0 0 M03
0 0 0 0
0 0 0 0
M03 0 0 0
//Migration matrix 3
0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
//Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
3 historical events
TcomP 2 3 1 1 0 1
TcomD 1 0 1 1 0 2
Tpyra 3 0 1 1 0 3
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 3.9e-8 OUTEXP
