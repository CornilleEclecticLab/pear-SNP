#!/usr/bin/env bash

#SBATCH -J s02.zip_fastStructure_and_run_CLUMPAK.sh
#SBATCH -o s02.zip_fastStructure_and_run_CLUMPAK.NOGIT.%J.out
#SBATCH -e s02.zip_fastStructure_and_run_CLUMPAK.NOGIT.%J.err
#SBATCH -c 1
#SBATCH --mem=2G

# @File     :   s02.zip_fastStructure_and_run_CLUMPAK.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/01/08 10:21:20
#Description:


####################
# load config file
####################
source s00.config.sh


####################
# zip FastStructure output
####################
random_id=${RANDOM}

cd ${OUTPUT_DIR}

INPUT_CLUMPAK_ZIP="input_clumpak_${random_id}.zip"

if [ -f $INPUT_CLUMPAK_ZIP ]; then
    rm $INPUT_CLUMPAK_ZIP
fi

mkdir result && cd result

cp ../K*/*.meanQ .

for i in `seq $MIN_K $MAX_K`
    do zip -q K${i}.zip *.${i}.meanQ
    done

zip -q ../${INPUT_CLUMPAK_ZIP} K*.zip

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

# This is unnecessary for on IFB cluster
# export LC_ALL=C


if [ -d "tmp_clumpak_${random_id}" ]; then
    rm -r "tmp_clumpak_${random_id}"
fi

OUTPUT_CLUMPAK="./output_clumpak_${random_id}"

if [ -d "${OUTPUT_CLUMPAK}" ]; then
    rm -r "${OUTPUT_CLUMPAK}"
fi

# NOTE: "./" in the --dir paramater is necessary
./CLUMPAK.pl \
    --id "tmp_clumpak_${random_id}" \
    --dir "${OUTPUT_CLUMPAK}" \
    --file "${OUTPUT_DIR}/${INPUT_CLUMPAK_ZIP}" \
    --inputtype admixture \
&& mv "${OUTPUT_CLUMPAK}" "${OUTPUT_DIR}/" 
