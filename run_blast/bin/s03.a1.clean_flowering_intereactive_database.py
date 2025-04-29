#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s03.a1.clean_flowering_intereactive_database.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/04/28 18:30:23
# @Description:
#    

import datetime
import pandas as pd
import textwrap
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')

tables = pd.read_html('../database/Flowering Interactive Database [FLOR-ID] - Flowering time.html')

df = tables[0]

df.to_csv("../database/Flowering_Interactive_Database.flowering_data.tsv", index=False)


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))