
cd /work/zruilin/pear/run_FastStructure/output/batch2.removeCloneAnd3LowQualEuropeWildPear
mkdir result
cp ./K*/*.meanQ ./result

cd result
for i in `seq 2 15`
do
    zip -q K${i}.zip *.${i}.meanQ
done

zip -q ../input_clumpak.zip K*.zip

rm -rf K*.zip
