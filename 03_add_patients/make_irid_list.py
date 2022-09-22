import numpy as np
import pandas as pd

dn='/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing'

list1 = pd.read_csv(f'{dn}/list1_bcdiag.csv', index_col=0)
list2 = pd.read_csv(f'{dn}/list2_bcdiag_provenc.csv', index_col=0)
list3 = pd.read_csv(f'{dn}/list3_bcdiag_provenc_gennote.csv', index_col=0)

ir_ids = []

for patlist in [list1, list2, list3]:
    ir_ids = ir_ids + list(patlist['ir_id'].values)

with open(f'{dn}/irid_list_for_sql_query.txt', "w") as f:
    for i in range(len(ir_ids)):
        if (i+1) % 10 == 0:
            f.write("(" + str(ir_ids[i]) + ")\n")
        else:
            f.write("(" + str(ir_ids[i]) + "), ")

