#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

sens_file = snakemake.input.sens_file
anti_file = snakemake.input.anti_file
out_file = str(snakemake.output)

sens_data = pd.read_csv(sens_file, sep='\t', comment='#', index_col=0)
anti_data = pd.read_csv(anti_file, sep='\t', comment='#', index_col=0)
data = pd.DataFrame(columns=sens_data.columns, index=sens_data.index)
data.index = sens_data.index

col_info = ['Chr', 'Start', 'End', 'Strand', 'Length']

sens = 0
anti = 0
for col in sens_data.columns:
    if col in col_info:
        data[col] = sens_data[col]
    else: 
        print(f'{col}')
        print(f'Sens: {np.sum(sens_data.loc[:, col])}')
        print(f'Anti: {np.sum(anti_data.loc[:, col])}')
        sens += np.sum(sens_data.loc[:, col])
        anti += np.sum(anti_data.loc[:, col])
        data[col] = [sens_data.loc[i, col] - anti_data.loc[i, col] for i in data.index ]
        print(f'Minimum: {min(data[col])}')
        print(f'Maximum: {max(data[col])}')
        data[col][data[col] <  0] = 0 
        print('------------------')

print(f'Sens: {sens}')
print(f'Anti: {anti}')
data.to_csv(out_file, sep='\t')
