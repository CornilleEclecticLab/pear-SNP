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

OUTPUT_CLUMPAK="./output_clumpak"

if [ -d "${OUTPUT_CLUMPAK}" ]; then
    rm -r "${OUTPUT_CLUMPAK}"
fi

#Xilong: If you are using the anther language, the CLUMPAK may have some bug, LC_ALL=C is changing the language to C(C for computer)
export LC_ALL=C

random_id=${RANDOM}

if [ -d "tmp_clumpak_${random_id}" ]; then
    rm -r "tmp_clumpak_${random_id}"
fi

# NOTE: "./" in the --dir paramater is necessary
./CLUMPAK.pl \
    --id "tmp_clumpak_${random_id}" \
    --dir "${OUTPUT_CLUMPAK}" \
    --file "${OUTPUT_DIR}/input_clumpak.zip" \
    --inputtype admixture \
&& \
mv "${OUTPUT_CLUMPAK}" "${OUTPUT_DIR}/"
