# Run_CLUMPAK_on_NEOcluster
Software CLUMPAK generally runs on its websites, but it is slow and limited by the network environment. CLUMPAK also provides a local version, I installed it in NEO and you can run it locally.

Since this software is a combination of a large number of perl scripts, there are a lot of dependencies between these scripts. I modified the main program's script to add some dependency libraries. But it can only run in the directory where the software is located. I changed the folder permissions and now everyone can read and write files under this folder.

Here are the pipeline steps:

```shell
cd /home/dygap/xchen/software/CLUMPAK/26_03_2015_CLUMPAK/CLUMPAK
export LC_ALL=C # If you are using the anther language, the CLUMPAK may have some bug, LC_ALL=C is changing the language to C(C for computer)
perl CLUMPAK.pl --id <job_name> --dir <output_dir> --file <your CLUMPAK_zip file>
## if your input file are genrate from admixture or fastStructure
perl CLUMPAK.pl --id <job_name> --dir <output_dir> --file <your CLUMPAK_zip file> --inputtype admixture
```
You can give any name for <job_name>, eg “mydata123”. But the <output_dir> must be a folder in CLUMPAK folder, eg ./output_test
But the output_dir should be different from job_name

The job will run about 3 to 4 hours on a K2 to K15 by repeating 20 runs of a dataset containing thirteen markers among 800 individuals. 
When the jobs are finished, you can copy the <output_dir> to your own folder. And please **remove** the <output_dir> locally.

PS, This pipeline is the “Main pipeline” of CLUMPAK, it also have Distruct for many K's, Compare and Best K pipeline. You can find more details on its website.

Reference
http://clumpak.tau.ac.il/
"CLUMPAK: a program for identifying clustering modes and packaging population structure inferences across K". Kopelman, Naama M; Mayzel, Jonathan; Jakobsson, Mattias; Rosenberg, Noah A; Mayrose, Itay. Molecular Ecology Resources 15(5): 1179-1191, doi: 10.1111/1755-0998.12387
http://clumpak.tau.ac.il/download/CLUMPAK_Documentation.pdf


st2ad.pl
Update 2022
if the STRUCTURE result didn't have Lable column, the CLUMPAK will fail
This script can convert it to admixture result
the output file will be name to input_file.Q
