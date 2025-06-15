import sys
import time
from collections import defaultdict
from collections import Counter

import pandas as pd


PRIORITY = """\
    miRNA rRNA tRNA Mt-tRNA snoRNA scRNA srpRNA snRNA antisense lncRNA RNA misc_RNA
    Cis-reg ribozyme RC IRES frameshift_element
    LINE SINE Simple_repeat Low_complexity Satellite DNA LTR CDS UTR5 UTR3 ncRNA
    intron Other Unknown antisense_rRNA\
""".split()
TYPE_TO_NUM = defaultdict(lambda: PRIORITY.index('Other'))
for i, rType in enumerate(PRIORITY):
    TYPE_TO_NUM[rType] = i
NUM_TO_TYPE = {v: k for k, v in TYPE_TO_NUM.items()}
PRIORITY_LEN = len(PRIORITY)

# def pri(annoSet):
#     for rType in PRIORITY:
#         if rType in annoSet: 
#             return rType
#     return 'Other'

# def file_open(fileIn):
#     if fileIn.split('.')[-1] == 'gz':
#         with gzip.open(fileIn, 'rt') as i:
#             file = i.readlines()
#         print("file_open : complete")
#         return file
#     if fileIn.split('.')[-1] != 'gz':
#         with open(fileIn) as i:
#             file = i.readlines()
#         print("file_open : complete")
#         return file

def main(args):
    fileIn = args[1]
    fileOut = args[2]
    # file = file_open(fileIn)
    # file_len = len(file)
    dReadAnnos = defaultdict(lambda: PRIORITY_LEN)
    # iterCount = 0
    # lstPercent = list(range(0, 101, 10))

    if fileIn.endswith('.gz'):
        import gzip
        file_open = gzip.open
    else:
        file_open = open

    with file_open(fileIn, 'rt') as fIn:
        for line in fIn:
            lstField = line.strip().split('\t')
            readName = lstField[3].split('/')[0]
            readAnno = lstField[-3].split('|')[0].replace('?', '')
            priority = TYPE_TO_NUM[readAnno]
            if dReadAnnos[readName] > priority:
                dReadAnnos[readName] = priority
            # dReadAnnos[readName].add(readAnno)
            # iterCount += 1
            # if int(iterCount/file_len * 100) in lstPercent: 
            #     print('parsing progress : ', int(iterCount/file_len * 100), '%')
            #     lstPercent.pop(0)
    print(time.ctime(), f" --- line parsing : complete --- {fileIn}")
    dctCounts = Counter(dReadAnnos.values())
    dctCounts = {NUM_TO_TYPE[k]: v for k, v in dctCounts.items()}
    seAnnotype = pd.Series(dctCounts)
    print(time.ctime(), f" --- series generation : complete --- {fileIn}")
    seAnnotype.to_csv(fileOut, sep='\t', header=None)

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv)
