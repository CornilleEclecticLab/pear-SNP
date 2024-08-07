#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s01.x1.format_mdist.py
# @Version  : 1.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/01/12 17:29:54
# @Description:
#     


with open("plink.mdist", "r") as fm:
    with open("plink.mdist.id","r") as fi:
        with open("plink.formatted.mdist", "w") as fo:
            len_id = len(fi.readlines())
            fi.seek(0)
            fo.write(f"\t{len_id}\n")
            for lm,li in zip(fm,fi):
                fo.write(f"{li.strip().split()[-1]}\t{lm.strip()}\n")
