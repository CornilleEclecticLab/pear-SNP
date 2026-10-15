//Parameters  for the coalescence simulation program : simcoal.exe
4 samples to simulate :
//Population effective sizes (number of genes)
Nbetu
Npash
NSand_CN
NWhite
//Samples sizes, samples ages, and inbreeding levels
16
36
52
12
//Growth rates  : negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
1
//Migration matrix 0
0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
//Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
3 historical events
TWhite 3 2 1 1 0 0
TSand_CN 2 1 1 1 0 0
Tbetu 0 1 1 1 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 3.9e-8 OUTEXP
