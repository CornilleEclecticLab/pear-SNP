#!/usr/bin/env bash

# @File     :   s03.CV_plot.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/03/12 16:14:54
#Description:
    # 


source s00.config.sh

${LOAD_PERL5}
${LOAD_R}

cd $OUTPUT_DIR

mkdir -p log

cp K*/*.log log

perl "${WORK_DIR}/bin/s03.x1.fastst_cv.pl"

Rscript "${WORK_DIR}/bin/s03.x2.cross-entropy.fig.R"

rm -r log
