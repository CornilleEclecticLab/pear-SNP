
cd /work/zruilin/Peach/run_FastStructure/output/s01.generate_fast
for rep in `seq 2 10`
do
        cd ./K$rep
        ls ./ | xargs -I {} sbatch -c 1 --mem=1G {}
        cd ..
done