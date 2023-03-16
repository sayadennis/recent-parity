import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import forestplot as fp

import string # only for troubleshooting purposes 

din = '/share/fsmresfiles/breast_cancer_pregnancy/stat_results'
dout = '/share/fsmresfiles/breast_cancer_pregnancy/plots'

results_fn = [
    # mutation vs parity 
    '20230307_mutation_vs_recencyparity5.csv', 
    '20230307_mutation_vs_recencyparity10.csv', 
    '20230307_mutation_vs_parity.csv',
    # tumor characteristics vs parity 
    '20230307_tumchar_vs_recencyparity5.csv', 
    '20230307_tumchar_vs_recencyparity10.csv', 
    '20230307_tumchar_vs_parity.csv',
]

for fn in results_fn:
    df = pd.read_csv(f'{din}/{fn}')

    # df = pd.DataFrame(
    #     [['all', 1.4, 0.9, 1.7],
    #      ['BRCA1', 1.2, 0.3, 1.4],
    #      ['BRCA2', 1.8, 1.2, 3.3]],
    #     index=['all', 'BRCA1', 'BRCA2'],
    #     columns=['var', 'or', 'low', 'high']
    # )

    # #### troubleshooting purposes only ####
    df = df.iloc[~pd.isnull(df['or']).values,:]
    # df['varname'] = list(df.index) # list(string.ascii_lowercase)[:df.shape[0]]
    # #######################################

    ax = fp.forestplot(df,  # the dataframe with results data
        estimate='or',  # col containing estimated effect size 
        ll='low', hl='high',  # columns containing conf. int. lower and higher limits
        varlabel='varname',  # column containing variable label
        ylabel='Confidence interval',  # y-label title
        xlabel='Odds Ratio',  # x-label title
        color_alt_rows=True,
        figsize=(4,4),
        **{
            'xline' : 1.,
            # "marker": "D",  # set maker symbol as diamond
            # "markersize": 35,  # adjust marker size
            "xlinestyle": (0, (10, 5)),  # long dash for x-reference line 
            # "xlinecolor": "#808080",  # gray color for x-reference line
            # "xtick_size": 12,  # adjust x-ticker fontsize
        }
    )
    ax.set_ylim(-0.5, df.shape[0]-0.5)

    # datestring = datetime.now().strftime("%Y%m%d")
    fout = 'forest_plot_' + fn.split('.')[0] + '.png'
    plt.savefig(f'{dout}/{fout}', bbox_inches='tight')
    plt.close()
