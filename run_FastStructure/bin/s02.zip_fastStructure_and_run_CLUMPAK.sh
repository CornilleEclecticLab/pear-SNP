#!/usr/bin/env bash

#SBATCH -J s02.zip_fastStructure_and_run_CLUMPAK.sh
#SBATCH -o s02.zip_fastStructure_and_run_CLUMPAK.NOGIT.%J.out
#SBATCH -e s02.zip_fastStructure_and_run_CLUMPAK.NOGIT.%J.err
#SBATCH -c 1
#SBATCH --mem=16G

# @File     :   s02.zip_fastStructure_and_run_CLUMPAK.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/01/08 10:21:20
#Description:



# load config file
source s00.config.sh


####################
# zip FastStructure output
####################
cd ${OUTPUT_DIR}

mkdir result && cd result

cp ../K*/*.meanQ .

for i in `seq $MIN_K $MAX_K`
    do zip -q K${i}.zip *.${i}.meanQ
    done

zip -q ../input_clumpak.zip K*.zip

rm K*.zip
rm *.meanQ

cd ${OUTPUT_DIR}
rmdir ./result


####################
# run CLUMPAK
####################
# load perl5
${LOAD_PERL5}

cd ${CLUMPAK_DIR}

random_id=${RANDOM}
# NOTE: "./" in the --dir paramater is necessary
./CLUMPAK.pl \
    --id "tmp_clumpak_${random_id}" \
    --dir ./output_clumpak \
    --file "${OUTPUT_DIR}/input_clumpak.zip" \
    --inputtype admixture

mv "./output_clumpak" \
    "${OUTPUT_DIR}/"
