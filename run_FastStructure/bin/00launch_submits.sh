
cd /work/zruilin/pear/run_FastStructure/output/batch3.removeLoquatColoneAnd3LowQualEuropeWildPear
for rep in `seq 2 15`
do
        cd ./K$rep
        ls ./ | xargs -I {} sbatch -c 1 --mem=1G {}
        cd ..
done