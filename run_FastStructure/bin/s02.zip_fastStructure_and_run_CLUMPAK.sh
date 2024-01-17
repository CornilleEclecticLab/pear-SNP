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
# Load config file
####################
source s00.config.sh


####################
# Zip FastStructure output
####################
cd "${OUTPUT_DIR}" || exit 1

random_id=${RANDOM}
INPUT_CLUMPAK_ZIP="input_clumpak_${random_id}.zip"

while [ -f "$INPUT_CLUMPAK_ZIP" ]; do
    random_id=${RANDOM}
    INPUT_CLUMPAK_ZIP="input_clumpak_${random_id}.zip"
done

TEMP_RESULT_DIR="temp_result_${random_id}"

if [ -d "$TEMP_RESULT_DIR" ]; then
    rm -r "$TEMP_RESULT_DIR" || exit 1 # exit with error code 1
fi

mkdir "${TEMP_RESULT_DIR}" && cd "${TEMP_RESULT_DIR}"

cp ../K*/*.meanQ .

for i in $(seq "$MIN_K" "$MAX_K"); do
    zip -q "K${i}.zip" *.${i}.meanQ
done

zip -q "../${INPUT_CLUMPAK_ZIP}" K*.zip

rm K*.zip
rm *.meanQ

cd ${OUTPUT_DIR} || exit 1
rmdir "${TEMP_RESULT_DIR}"


####################
# run CLUMPAK
####################
# Load perl5 and ghostscript
${LOAD_PERL5}
${LOAD_GHOSTSCRIPT}

cd "${CLUMPAK_DIR}" || exit 1

# This is unnecessary for on IFB cluster
# export LC_ALL=C

CLUMPAK_JOB="clumpak_${random_id}"
if [ -d "$CLUMPAK_JOB" ]; then
    rm -r "$CLUMPAK_JOB" || exit 1
fi

OUTPUT_CLUMPAK="./output_clumpak_${random_id}"
if [ -d "${OUTPUT_CLUMPAK}" ]; then
    rm -r "${OUTPUT_CLUMPAK}" || exit 1
fi

# NOTE: "./" in the --dir paramater is necessary
./CLUMPAK.pl \
    --id "${CLUMPAK_JOB}" \
    --dir "${OUTPUT_CLUMPAK}" \
    --file "${OUTPUT_DIR}/${INPUT_CLUMPAK_ZIP}" \
    --inputtype admixture \
&& mv "${OUTPUT_CLUMPAK}" "${OUTPUT_DIR}/" 
