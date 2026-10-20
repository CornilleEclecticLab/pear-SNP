#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     :   find_admixed_individuals.py
# @Version  :   2.1.0
# @Author   :   NIE Yuqi
# @Time(CET):   2023-04-14 13:49:26
#Description:
    # Update: 2024-02-15 16:31:12
    # The color list is updated, and specific colors for different dataset.
    
    # Update: v2.1.0
    # add color palette for European mapped onto P.communis

import datetime
import sys

#################  Parameters  #################
threshold = 0.8 
print(f'Threshold for pure samples is set as: {threshold}') 


## Paste your colors here, and consistent with the order of you used in pophelper
# col_F = ("#648fff",
#          "#785ef0",
#          "#dc267f",
#          "#fe6100",
#          "#ffb000",
#          "#FFA07A",  # (Light Salmon)
#          "#FF69B4",  # (Hot Pink)
#          "#1E90c8",  # changed from #1E90ff
#          "#DC143C",  # (Crimson)
#          "#8B008B",  # (Dark Magenta)
#          "#32CD32",  # (Lime Green)
#          "#00CED1",  # (Dark Turquoise)
#          "#FF1493",  # (Deep Pink)
#          "#00FF7F",  # (Spring Green)
#          "#9ACD32",  # (Yellow Green)
#          "#FF00FF",  # (Magenta)
#          "#4B0082",  # (Indigo)
#          "#FFD700"  # (Gold)
# ) # 2024-02-15-Full

# col_A = ("#648fff",  # Full-1
#          "#dc267f",  # Full-3
#          "#fe6100",  # Full-4
#          "#ffb000",  # Full-5
#          "#FF69B4",  # Full-7
#          "#DC143C",  # Full-9
#          "#8B008B",  # Full-10
#          "#32CD32",  # Full-11
#          "#00CED1",  # Full-12
#          "#FF1493",  # Full-13
#          "#00FF7F",  # Full-14
#          "#9ACD32",  # Full-15
#          "#FF00FF",  # Full-16
#          "#4B0082",  # Full-17
#          "#785ef0",  # Full-2
#          "#FFA07A",  # Full-6
#          "#1E90c8",  # Full-8
#          "#FFD700"  # Full-18
# ) # 2024-02-15-ASI

# col_E = ("#785ef0",  # Full-2
#          "#FFA07A",  # Full-6
#          "#1E90c8",  # Full-8
#          "#FFD700",  # Full-18
#          "#8B008B",  # Full-10
#          "#32CD32",  # Full-11
#          "#00CED1",  # Full-12
#          "#FF1493",  # Full-13
#          "#00FF7F",  # Full-14
#          "#9ACD32",  # Full-15
#          "#FF00FF",  # Full-16
#          "#4B0082",  # Full-17
#          "#648fff",  # Full-1
#          "#dc267f",  # Full-3
#          "#fe6100",  # Full-4
#          "#ffb000",  # Full-5
#          "#FF69B4",  # Full-7
#          "#DC143C"  # Full-9
# ) # 2024-02-15-EUO


# col_F = [
#     "#960558",
#     "#442CF4",
#     "#FFA07A",
#     "#5E4734",
#     "#FF5800",
#     "#27bc40",
#     "#FFB727",
#     "#E599F7",
#     "#9A51FF",
#     "#00CED1",
#     "#648FFF",
#     "#1E90c8",
#     "#FF69B4",
#     "#DC143C",
#     "#8B008B",
#     "#FF00FF",
#     "#80a492",
#     "#035b03"
# ]

# col_A = [
#     "#960558",
#     "#FFA07A",
#     "#5E4734",
#     "#FF5800",
#     "#FFB727",
#     "#FFC9C9",
#     "#FD79A8",
#     "#B2F2BB",
#     "#91FBEF",
#     "#99daf7",
#     "#27bc40",
#     "#FF69B4",
#     "#4D2547",
#     "#035b03",
#     "#DC143C",
#     "#648FFF"
# ]

# col_E = [
#     "#442CF4",
#     "#27bc40",
#     "#E599F7",
#     "#035b03",
#     "#9b59b6",
#     "#DC143C",
#     "#648FFF",
#     "#00CED1",
#     "#1E90c8",
#     "#8B008B",
#     "#9A51FF",
#     "#960558",
#     "#FFA07A",
#     "#5E4734",
#     "#FF5800",
#     "#FFB727"
# ]
# color_pallet_version = '# 2024-02-28'

col_F = ["#960558",
         "#442CF4",
         "#FFA07A",
         "#5E4734",
         "#FF5800",
         "#27bc40",
         "#FFC9C9",
         "#91FBEF",
         "#FFB727",
         "#9A51FF",
         "#FD79A8",
         "#B2F2BB",
         "#E599F7",
         "#99daf7",
         "#00CED1",
         "#648FFF",
         "#1E90c8",
         "#FF69B4",
         "#DC143C",
         "#8B008B",
         "#FF00FF",
         "#80a492",
         "#035b03" 
]

col_A = ["#960558",
         "#FFA07A",
         "#5E4734",
         "#FF5800",
         "#FFB727",
         "#FFC9C9",
         "#FD79A8",
         "#E599F7",
         "#99daf7",
         "#9A51FF",
         "#27bc40",
         "#FF69B4",
         "#4D2547",
         "#035b03",
         "#DC143C",
         "#648FFF" 
]

col_E = ["#442CF4",
         "#27bc40",
         "#91FBEF",
         "#035b03",
         "#b2f2bb",
         "#DC143C",
         "#648FFF",
         "#00CED1",
         "#1E90c8",
         "#8B008B",
         "#9b59b6",
         "#960558",
         "#FFA07A",
         "#5E4734",
         "#FF5800",
         "#FFB727"
]


col_EE = ["#442CF4",  # Full-1  # blue
         "#27bc40",  # Full-5  # green
         "#91FBEF",  # Full-9 
         "#035b03",  # Full-15 # dark green 2024-02-28 11:35:26 changed to a bit lighter color from #033303
         "#648FFF",  # Full-7  # sky blue
         "#b2f2bb",
         "#DC143C",  # Full-12
         "#8B008B",  # Full-11
         "#9b59b6",  # Full-14
         "#960558",  # Full-0  # red
         "#FFA07A",  # Full-2  # vived orange
         "#1E90c8",  # Full-8 
         "#00CED1",  # Full-13
         "#5E4734",  # Full-3  # brown
         "#FF5800",  # Full-4  # orange
         "#FFB727"  # Full-6  # yellow
 ] # 2024-08-04-EUO-For_mapping_on_comm
   #  From branch17.batch45

color_pallet_version = '# 2025-8-04'

try:
    if sys.argv[1] != '' and sys.argv[2] != '':
        pass
except IndexError:
    print(f'Usage: python3 {sys.argv[0]} <input_Qfile> <F|A|E|EE>')
    sys.exit(1)

start_time = datetime.datetime.now()
print("{0:=^40}".format(' Start '))

input_file = sys.argv[1]
fo = open(input_file+"admixed.tsv", 'w', encoding='utf-8')
fo2 = open(input_file+"admixed_color.tsv", 'w', encoding='utf-8')

if sys.argv[2] in list('FAE') or sys.argv[2] == 'EE':
    color_ls = eval("col_"+sys.argv[2])
else:
    print('The second parameter expected to be F, A or E')
    sys.exit(1)

with open(input_file, 'r') as fi:
    n = 0
    count_pure = 0
    for line in fi:
        admix = True
        n+=1
        line = line.strip()
        Qs = line.split()
        for i, value in enumerate(Qs):
            if float(value) >= threshold:
                count_pure += 1
                admix = False
                color = color_ls[i]
                break

        if admix:
            color = '#888888'
            out = "Sample"+str(n)+"\t"+"\t".join(Qs)
            print(out)
            fo.write(out+'\n')

        fo2.write(color+'\n')

print(f'Number of pure samples: {count_pure}')

end_time = datetime.datetime.now()
print('')
print(' END '.center(40,'='))
print(str(end_time-start_time).center(40))
