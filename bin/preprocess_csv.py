'''

This script is used to preprocess the data. 
The fixation fraction calculated online is replaced with the fixation fraction calculated offline,
for the experiment 3.
The data on particular date and with fixation fraction < 89 and removed.
###
'''
import numpy as np
import pandas as pd
# from pathlib import Path


new_names = {'background': 'BG',
             'ori_faces': 'DF',
             'ori_bodies': 'DB',
             'ori_objects': 'DO',
             'scram_faces': 'MF',
             'scram_bodies': 'MB',
             'scram_objects': 'MO',
             'pscram_faces': 'PF',
             'pscram_bodies': 'PB',
             'pscram_objects': 'PO',
             'stat_faces': 'SF',
             'stat_bodies': 'SB',
             'stat_objects': 'SO', }

expt_new_names = {'Dyn_mScram': 'Expt. 1',
                  'Dyn_pScram': 'Expt. 2',
                  'Dyn_stat': 'Expt. 3'
                  }


# load data from each date, monkey, and expt


def concat_data(data_path):
    sc_paths = data_path.glob(f'**/data_table*.csv')
    sc_paths = [x for x in sc_paths]
    sc_paths.sort()  # sort by date, monkey, and expt

    # concat all data
    df_all = []
    for path in sc_paths:
        df = pd.read_csv(path)
        col_order = ['exp_name', 'monkey', 'date', 'run', 'condt', 'xfixwin', 'yfixwin', 'fixPerc', 'fixPerc_condt', 'srv_perc', 'xdata', 'ydata',
                     'Xvel', 'Yvel', 'xvel_flt', 'yvel_flt', 'sacc', 'dist']
        df = df.reindex(columns=col_order)
        df['condt'] = df['condt'].map(new_names)
        df['exp_name'] = df['exp_name'].map(expt_new_names)
        df_all.append(df)
    df_all = pd.concat(df_all)
    df_all.reset_index(inplace=True)
    return df_all


def removeInvalidData(df_all):
    # replacing online fixtation fraction value with the estimated offline value in the experiment 3
    # estimat and remove data with fixation fraction calulted offline

    # get the fixation data from the eye movement
    df_st_smp = df_all[['xdata', 'ydata', 'exp_name',
                        'monkey', 'date', 'condt', 'run']]
    thrs_val = 1.5  # threshold value for 3 degree of visual angle
    df_st_cl = df_st_smp[(df_st_smp.xdata.abs() < thrs_val) &
                         (df_st_smp.ydata.abs() < thrs_val)]

    df_fix = df_st_cl.groupby(
        ['monkey', 'date', 'condt', 'run', 'exp_name']).count().reset_index()
    # df_fix.rename(columns={'xdata': 'fixPerc_new'}, inplace=True)
    df_ori = df_st_smp.groupby(
        ['monkey', 'date', 'condt', 'run', 'exp_name']).count().reset_index()

    # new fixation data
    df_fix['fix_new'] = np.divide(df_fix.xdata.values, df_ori.xdata.values)*100
    df_final = df_fix.groupby(
        ['exp_name', 'monkey', 'date', 'condt', 'run']).mean().reset_index()
    df_final.drop(columns=['xdata', 'ydata'], inplace=True)

    # merge df_all and df_final
    df_all = df_all.merge(
        df_final, on=['exp_name', 'monkey', 'date', 'condt', 'run'], how='left')

    # find the date and run with fixation fraction < 89 in expriment 3
    df_stt = df_final.query('exp_name == "Expt. 3"')
    df_stt_all = df_stt.groupby(
        ['date', 'exp_name', 'monkey', 'run']).mean().reset_index()
    df_rmn = df_stt_all.copy()
    df_rmn = df_rmn[df_rmn.fix_new.round() < 89]
    df_rmn['date_run'] = df_rmn.date.apply(str) + '_' + df_rmn.run.apply(str)

    # remove the data on the particular date
    df_raDate = df_all[df_all.date == 211019]
    df_all.drop(df_raDate.index, inplace=True)

    # remove the data on the particular date and run
    df_all['date_run'] = df_all.date.apply(str) + '_' + df_all.run.apply(str)
    df_rm = df_all[(df_all.date_run.isin(df_rmn.date_run))]
    df_all.drop(df_rm.index, inplace=True)

    # replace the fixation fraction with the new one for the experiment 3
    df_all.loc[df_all['exp_name'] == 'Expt. 3',
               'fixPerc_condt'] = df_all.loc[df_all['exp_name'] == 'Expt. 3', 'fix_new']
    df_all.to_csv('../results/data/df_all.csv', index=False)
    return df_all
