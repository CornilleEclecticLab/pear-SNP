
cd /work/zruilin/Peach/run_FastStructure/output/s01.generate_fast
mkdir result
cp ./K*/*.meanQ ./result

cd result
for i in `seq 3 4`
do
    zip -q K${i}.zip *.${i}.meanQ
done

zip -q ../input_clumpak.zip K*.zip

rm -rf K*.zip
