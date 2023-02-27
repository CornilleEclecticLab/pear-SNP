
cd /work/zruilin/pear/run_FastStructure/output/batch6.only_Prifolia
for rep in `seq 2 15`
do
        cd ./K$rep
        ls ./ | xargs -I {} sbatch -c 1 --mem=1G {}
        cd ..
done