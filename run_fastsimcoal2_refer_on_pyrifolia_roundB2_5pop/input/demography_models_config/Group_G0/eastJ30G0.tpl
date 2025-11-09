//Parameters  for the coalescence simulation program : simcoal.exe
5 samples to simulate :
//Population effective sizes (number of genes)
Nbetu
Npash
Npyri_JP
NSand_CN
NSandCNSE
//Samples sizes, samples ages, and inbreeding levels
16
34
36
54
20
//Growth rates  : negative growth implies population expansion
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
1
//Migration matrix 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
4 historical events
T2 2 1 1 1 0 0
T4 4 3 1 1 0 0
T3 3 1 1 1 0 0
T0 0 1 1 1 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 3.9e-8 OUTEXP
