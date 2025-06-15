from collections import defaultdict
import os
import sys

import pandas as pd

def main(args):
    stats_dir = args[1]
    export_file_path = args[2]
    file_list = os.listdir(stats_dir)
    stats_dict = parse_stat_files(stats_dir, file_list)
    pd.DataFrame(stats_dict).to_csv(export_file_path, sep='\t')

def parse_stat_files(stats_dir, file_list):
    stats_dict = defaultdict(dict)
    for file in file_list:
        items = file.split('.')
        assert items[2] == 'count'
        assert items[3] == 'txt'
        sample, category = items[0], items[1]
        with open(os.path.join(stats_dir, file), 'rt') as fIn:
            stats_list = fIn.readlines()
        assert stats_list
        if len(stats_list) == 1:
            stats_dict[sample][category] = int(stats_list[0])
        else:
            annotated = 0
            for stats in stats_list:
                rType, counts = stats.strip().split('\t')
                counts = int(counts)
                annotated += counts
                stats_dict[sample][rType] = counts
            stats_dict[sample]['annotated'] = annotated
    return stats_dict

if __name__ == '__main__':
    main(sys.argv)