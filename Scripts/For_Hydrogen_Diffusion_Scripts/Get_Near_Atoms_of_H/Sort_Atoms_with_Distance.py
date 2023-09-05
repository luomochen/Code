# 读取数据为DataFrame格式，对distance参数进行排序，并覆写原txt文件。
import pandas as pd
import sys

filename = 'Atom_distances/Atom_distances_0' + sys.argv[1] + '.txt'
Atom_info = pd.read_table(filename,sep='\s+')
Atom_distances_sort = Atom_info.sort_values(by='Distance')
Atom_distances_sort.to_csv(sep='\t', index=False)
Atom_info.sort_values(by='Distance').to_csv(filename,sep='\t', index=False)