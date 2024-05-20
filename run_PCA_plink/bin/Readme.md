# Steps 
1. Edite `s00.conifg.sh` for input path
2. Run `s01.IFB.get_pca_plink.sh` Slurm script, using `sbatch s01.IFB.get_pca_plink.sh`
3. Edite input files, `id_map.txt` and `populations_group.tsv` in the `../input/` folder
   Note that the color code should use double quotes, e.g. "#888888". to be readable by R
4. cd to the `../output/<sub_folder>` folder
5. Run `s02.format_pca_data_for_plot.py` python script located in the `bin` folder
6. Run `s03.plot_PCA.R` R script


# PVE
PVE, percentage variance explained

Many tutorials uses the following formula:

pve = eigenval/sum(eigenval)*100

But as plink only shows the first 20 eigenval by default, if using this formula in the case not all eigenval present, we would overestimate the PVE

See https://www.biostars.org/p/436945/ for more information.

If we let plink output all the PCs, we will get the true sum of eigenval.

Note that, if let plink outputs more than 20 PCs, the minus sign of PC values may display incorrectly. It happened when PC10, during my test.
