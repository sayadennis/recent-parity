import os
import sys
import numpy as np
import pandas as pd
import random
import re

din = '/share/fsmresfiles/breast_cancer_pregnancy/notes/genetic_counsel'
dindv = '/share/fsmresfiles/breast_cancer_pregnancy/notes/genetic_counsel/individual_txt'
ddev = '/share/fsmresfiles/breast_cancer_pregnancy/notes/genetic_counsel/dev_set'

query = pd.read_csv(os.path.join(din, 'notes_genetic_counsel.csv'), header=None)

for pat_irid in query[0].unique():
    textlist = query.iloc[[x==pat_irid for x in query[0].values],6].values
    for i in range(len(textlist)):
        with open(os.path.join(dindv, str(pat_irid) + '_{}.txt'.format(str(i))), 'w') as f:
            f.write(textlist[i])

dev_patids = random.sample(list(query[0].unique()), 50)

for pat_irid in dev_patids:
    textlist = query.iloc[[x==pat_irid for x in query[0].values],6].values
    for i in range(len(textlist)):
        with open(os.path.join(ddev, str(pat_irid) + '_{}.txt'.format(str(i))), 'w') as f:
            f.write(textlist[i])

