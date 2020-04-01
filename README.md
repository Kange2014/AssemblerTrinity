# AssemblerTrinity
This is one Torrent Suite Software plugin for sequencing data assembly from ion torrent platforms, like genestudio S5 serires. It provides one intuitive and graphic interface to users for easy-use.

The basic workflow for this plugin is:
1. use trinity to de novo assemble reads from ion torrent sequencing platform directly with/without genome-guided mode;

2. each sample’s reads are then aligned to its de novo assembly using tmap. Variant positions in each assembly were identified using tvc on the read alignments. The assembly was refined to represent the major allele (allele frequency >= 51%) at each variant site. This align-call-refine cycle is iterated twice, to minimize reference bias in the assembly.
