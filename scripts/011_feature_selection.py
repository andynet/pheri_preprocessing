#!/usr/bin/python3

import sys
import pandas as pd
from sklearn.feature_selection import VarianceThreshold

if len(sys.argv) != 4:
    print('Usage:', sys.argv[0], '<matrix.raw.tsv> <matrix.fs.tsv> <variance_threshold>')
    exit()

matrix = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
selection = VarianceThreshold(threshold=(float(sys.argv[3]) * (1 - float(sys.argv[3]))))
indices = selection.fit(matrix).get_support(indices=True)
selected_matrix = matrix.iloc[:, indices]

selected_matrix.to_csv(sys.argv[2], sep='\t')
