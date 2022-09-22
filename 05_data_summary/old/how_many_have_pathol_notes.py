import os
import numpy as np
import pandas as pd
import re
from datetime import datetime
from matplotlib_venn import venn2
from matplotlib import pyplot as plt

din='/share/fsmresfiles/breast_cancer_pregnancy/data/pathology'
dout='/share/fsmresfiles/breast_cancer_pregnancy/data_summary'

notepath = pd.read_csv(f'{din}/notes_pathology.csv', header=None)
oncpath = pd.read_csv(f'{din}/pathology_from_onc_tables.csv', header=None)
