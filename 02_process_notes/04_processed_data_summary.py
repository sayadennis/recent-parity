# the purpose of this script is to summarize the data collected from the SQL query and note parsing for the recent parity project.
import os
import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data'

joined = pd.read_csv(os.path.join(dn, 'recent_parity_combined_data.csv'))

#############################
#### assess missing rate ####
#############################

def missing_rate(table):
    misstab = pd.DataFrame(index=['%missing'], columns=table.columns)
    for colname in table.columns:
        misstab.loc['%missing', colname] = np.sum(pd.isnull(table[colname]))*100/len(table[colname])
    return misstab

misstab = missing_rate(joined)

#########################################################
#### Plot summary visualizations for genetic testing ####
#########################################################

genelist = ['brca', 'brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm']

def plot_overall(table, outdir):
    import matplotlib as mpl
    mpl.rcParams['font.size'] = 14.0
    has_positive = 0
    has_variant = 0
    negative = 0
    for i in table.index:
        if np.any(table.loc[i,:]=='positive'):
            has_positive += 1
        elif np.any(table.loc[i,:]=='variant'):
            has_variant += 1
        else:
            negative +=1
    # plot pie chart
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(12,12))
    labels = ['Has positive', 'No positive but has variant', 'No positive or variant']
    sizes = [has_positive, has_variant, negative]
    axs.pie(sizes, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90, colors=['violet', 'cornflowerblue', 'yellowgreen'])
    axs.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    axs.set_title('Genetic testing results', fontsize=20)
    axs.legend(loc='lower left')
    fig.savefig(os.path.join(outdir, 'overall_genetic_testing_result_piechart.png'))
    return

def plot_gen(table, outdir):
    ## Plot the rate of "only generals" - i.e. only have "Positive" or "Negative" results with no gene names 
    only_general = np.sum(~pd.isnull(table['general']) & [np.all(pd.isnull(table.iloc[i,2:])) for i in range(table.shape[0])])
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    labels = ['General positive/negative only', 'Gene names available']
    sizes = [only_general, table.shape[0]-only_general]
    fig, ax = plt.subplots(figsize=(10,6))
    ax.pie(sizes, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90, colors=['orange', 'cornflowerblue'])
    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    fig.savefig(os.path.join(outdir, 'general_posneg_piechart.png'))
    ## Plot the positive/negative/missing for genetic testing results 
    import matplotlib as mpl
    mpl.rcParams['font.size'] = 14.0
    fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(24,12))
    genelist = ['brca', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2','atm']
    for i in range(8):
        if i==0:
            only_brca = np.sum(~pd.isnull(table['brca']) & [np.all(pd.isnull(table.iloc[i,3:5])) for i in range(table.shape[0])])
            brca12 = np.sum(
                np.array([table.loc[i,'brca1'] in ['positive', 'variant'] for i in range(table.shape[0])]) & 
                np.array([table.loc[i,'brca2'] in ['positive', 'variant'] for i in range(table.shape[0])])
            )
            brca1 = np.sum(
                np.array([table.loc[i,'brca1'] in ['positive', 'variant'] for i in range(table.shape[0])]) & 
                ~np.array([table.loc[i,'brca2'] in ['positive', 'variant'] for i in range(table.shape[0])])
            )
            brca2 = np.sum(
                ~np.array([table.loc[i,'brca1'] in ['positive', 'variant'] for i in range(table.shape[0])]) & 
                np.array([table.loc[i,'brca2'] in ['positive', 'variant'] for i in range(table.shape[0])])
            )
            neg = np.sum(
                np.array([table.loc[i,'brca']=='negative' for i in range(table.shape[0])]) &
                np.array([table.loc[i,'brca1']=='negative' for i in range(table.shape[0])]) &
                np.array([table.loc[i,'brca2']=='negative' for i in range(table.shape[0])])
            )
            miss = np.sum([np.all(pd.isnull(table.loc[i,['brca', 'brca1', 'brca2']])) for i in range(table.shape[0])])
            # plot
            labels = ['BRCA1&2', 'BRCA1', 'BRCA2', 'Negative', 'Missing'] # 'Positive for "BRCA"', 
            sizes = [brca12, brca1, brca2, neg, miss] # only_brca, 
            axs[0,i].pie(sizes, labels=labels, autopct='%1.1f%%',
                    shadow=True, startangle=90, colors=['lightcoral', 'orange', 'violet', 'cornflowerblue', 'yellowgreen'])
            axs[0,i].axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
            axs[0,i].set_title('BRCA')
            axs[0,i].legend(loc='lower left')
        else:
            varname = genelist[i]
            genename = varname.upper()
            pos = np.sum(np.array([table.loc[i,varname] in ['positive', 'variant'] for i in range(table.shape[0])]))
            neg = np.sum(np.array([table.loc[i,varname]=='negative' for i in range(table.shape[0])]))
            miss = table.shape[0] - pos - neg
            # plot
            labels = ['Positive/variant', 'Negative', 'Missing']
            sizes = [pos, neg, miss]
            if i<4:
                axs[0,i].pie(sizes, labels=labels, autopct='%1.1f%%',
                        shadow=True, startangle=90, colors=['violet', 'cornflowerblue', 'yellowgreen'])
                axs[0,i].axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
                axs[0,i].set_title(genename)
            else:
                axs[1,i-4].pie(sizes, labels=labels, autopct='%1.1f%%',
                        shadow=True, startangle=90, colors=['violet', 'cornflowerblue', 'yellowgreen'])
                axs[1,i-4].axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
                axs[1,i-4].set_title(genename)
    fig.savefig(os.path.join(outdir, 'genetic_test_piecharts.png'))
    return

plot_gen(joined, '/share/fsmresfiles/breast_cancer_pregnancy/data_summary')

#####################################
#### plot summary visualizations ####
#####################################

varnames = {
    'gravida' : 'Gravida',
    'para' : 'Para',
    'menarche' : 'Age of menarche',
    # 'menopause' : 'Manopausal status or Age',
    # 'lmp' : 'Last menstrual period',
    'agefirstpreg' : 'Age of first preg/birth',
    'ageoflastpreg' : 'Age of last preg/birth',
    'lastnursed' : 'Last nursed (years ago)',
    # 'duration' : 'Nursing duration'
}

def plot_parity_hists(table, outdir): # for parity data 
    misstab = missing_rate(table)
    for colname in table.columns:
        if colname in ['ir_id', 'title_type', 'lmp']:
            continue
        else:
            miss = misstab.loc['%missing', colname]
            data = table[colname]
            if colname in ['gravida', 'para', 'menarche', 'agefirstpreg', 'ageoflastpreg']: # histogram variables
                if colname in ['menarche']:
                    data = data.iloc[~pd.isnull(data.values) & np.array([x!='unknown' for x in data.values])].astype(float)
                fig, ax = plt.subplots()
                ax.hist(data.iloc[~pd.isnull(data).values])
                ax.set_xlabel(varnames[colname])
                ax.set_ylabel('Counts')
                if colname in ['menarche']:
                    ax.set_xlim([0,25])
                ax.set_title('Histogram of {} ({:.1f}% missing)'.format(varnames[colname], miss))
                fig.savefig(os.path.join(outdir, colname + '_histogram.png'))
                plt.close(fig)
            elif colname=='menopause': # bar graph + histogram variables
                fig, axs = plt.subplots(nrows=1, ncols=2)
                pre = np.sum(data=='premenopausal')
                postonly = np.sum(data=='postmenopausal')
                postage = data.iloc[
                    ~pd.isnull(data.values) 
                    & np.array([x!='premenopausal' for x in data.values])
                    & np.array([x!='postmenopausal' for x in data.values])
                    & np.array([x!='perimenopausal' for x in data.values])
                    & ~np.array([str(x).startswith('in') for x in data.values])
                ].astype(int)
                axs[0].bar(np.arange(3), [pre, postonly, len(postage)])
                axs[0].set_xticks(np.arange(3))
                axs[0].set_xticklabels(['pre', 'post', 'has_age'])
                axs[0].set_ylabel('Counts')
                axs[0].set_title('Menopausal status ({:.1f}% missing)'.format(miss))
                axs[1].hist(postage)
                axs[1].set_xlabel('Age of menopause')
                # axs[1].set_ylabel('Counts')
                axs[1].set_title('Menopausal age')
                fig.savefig(os.path.join(outdir, colname + '_histogram.png'))
                plt.close(fig)
            elif colname=='lastnursed':
                fig, axs = plt.subplots(nrows=1, ncols=2)
                yearsago = []
                for i in range(data.shape[0]):
                    try:
                        if int(data.iloc[i])<100:
                            yearsago.append(int(data.iloc[i]))
                    except:
                        continue
                fig, ax = plt.subplots()
                ax.hist(yearsago)
                ax.set_xlabel(varnames[colname])
                ax.set_ylabel('Counts')
                ax.set_title('Histogram of {} ({:.1f}% missing)'.format(varnames[colname], miss))
                fig.savefig(os.path.join(outdir, colname + '_histogram.png'))
                plt.close(fig)
    return

outdir = '/share/fsmresfiles/breast_cancer_pregnancy/data_summary'

plot_parity_hists(joined, outdir)

def plot_parity_missing(table, outdir):
    import matplotlib as mpl
    mpl.rcParams['font.size'] = 14.0
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20,12))
    plotlist = ['gravida', 'para', 'agefirstpreg', 'ageoflastpreg', 'lastnursed', 'has "lastpreg" or "lastnursed"']
    for i in range(6):
        if i<5: # everything other than the last plot
            varname = plotlist[i]
            hasvalue = np.sum(~pd.isnull(table[plotlist[i]]))
            miss = np.sum(pd.isnull(table[plotlist[i]]))
            labels = ['Has value', 'Missing']
            sizes = [hasvalue, miss]
            if i<3:
                axs[0,i].pie(sizes, labels=labels, autopct='%1.1f%%',
                    shadow=True, startangle=90, colors=['cornflowerblue', 'yellowgreen'])
                axs[0,i].axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
                axs[0,i].set_title(plotlist[i])
            else:
                axs[1,i-3].pie(sizes, labels=labels, autopct='%1.1f%%',
                    shadow=True, startangle=90, colors=['cornflowerblue', 'yellowgreen'])
                axs[1,i-3].axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
                axs[1,i-3].set_title(plotlist[i])
        else:
            varname = plotlist[i]
            hasvalue = np.sum(~pd.isnull(table['ageoflastpreg']) | ~pd.isnull(table['lastnursed']))
            miss = np.sum(pd.isnull(table['ageoflastpreg']) & pd.isnull(table['lastnursed']))
            axs[1,i-3].pie(sizes, labels=labels, autopct='%1.1f%%',
                shadow=True, startangle=90, colors=['cornflowerblue', 'yellowgreen'])
            axs[1,i-3].axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
            axs[1,i-3].set_title(plotlist[i])
    fig.savefig(os.path.join(outdir, 'parity_missing_piecharts.png'))
    return

plot_parity_missing(joined, '/share/fsmresfiles/breast_cancer_pregnancy/data_summary')
