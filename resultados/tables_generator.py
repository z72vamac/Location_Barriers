import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
# from scipy import stats
import matplotlib
import tikzplotlib
import os
import matplotlib.patches as mpatches
from scipy.stats import shapiro 
from scipy.stats import mannwhitneyu
import itertools

# Parametros de las representaciones graficas
sns.set(style="darkgrid")

matplotlib.rcParams["font.sans-serif"] = 'Times New Roman'
matplotlib.rcParams["font.serif"] = 'Times New Roman'
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['text.usetex'] = True


## Recogida de datos
pd.set_option('display.max_rows', 200)

results0 = pd.read_csv('resultados/results_0.csv').iloc[:, 10:]
results1 = pd.read_csv('resultados/results_1.csv').iloc[:, 10:]
results2 = pd.read_csv('resultados/results_2.csv').iloc[:, 10:]
results3 = pd.read_csv('resultados/results_3.csv').iloc[:, 10:]
results4 = pd.read_csv('resultados/results_4.csv').iloc[:, 10:]

parameters0 = pd.read_csv('resultados/parameters_0.csv').iloc[:, 2:]

results0 = pd.concat([parameters0, results0], axis=1)
results1 = pd.concat([parameters0, results1], axis=1)
results2 = pd.concat([parameters0, results2], axis=1)
results3 = pd.concat([parameters0, results3], axis=1)
results4 = pd.concat([parameters0, results4], axis=1)

results = pd.concat([results0, results1, results2, results3, results4]).reset_index()

results.rename(columns={'single':'approach'}, inplace=True)

results['perc_B'] = results['perc_B'].astype('category')
results['perc_B'] = results['perc_B'].cat.reorder_categories(['10.0%', '20.0%', '50.0%', '100%'])

results['k'] = results['k'].astype('category')
results['k'] = results['k'].cat.reorder_categories(['1', '10%', '25%'])

results['approach'] = results['approach'].astype('category')
results['approach'] = results['approach'].replace({True:'single', False:'multi'})
results['approach'] = results['approach'].cat.reorder_categories(['single', 'multi'])

results['init'] = results['init'].astype('category')
results['init'] = results['init'].replace({True:'with', False:'without'})
results['init'] = results['init'].cat.reorder_categories(['without', 'with'])

results['Gap'] *= 100
# results['ObjVal_h2'] = results[['ObjVal_h', 'ObjVal_h2']].min(axis=1)
# results.loc[results['init'] == True, 'Runtime'] = results.loc[results['init']== True, 'Runtime'].fillna(3600)
vector = ['Gap' , 'Runtime','NodeCount','ObjVal','Runtime_h','Runtime_h2','ObjVal_h','ObjVal_h2']

# results.loc[results['init'] == True, vector] = results.loc[results['init'] == True].groupby(['n_N', 'perc_B', 'k', 'approach', 'wL','lazy','A4']).transform(lambda x: x.fillna(x.mean())).reset_index().loc[:, vector]

results.loc[results['init'] == 'with', 'NodeCount'] = results.loc[results['init']== 'with', 'NodeCount'].ffill(limit=100)

results.loc[results['init'] == 'with', 'ObjVal_h'] = results.loc[results['init']== 'with', 'ObjVal_h'].ffill(limit=100)
results.loc[results['init'] == 'with', 'ObjVal_h2'] = results.loc[results['init']== 'with', ['ObjVal_h', 'ObjVal_h2']].min(axis=1)
results.loc[results['init'] == 'with', 'ObjVal'] = results.loc[results['init']== 'with', ['ObjVal', 'ObjVal_h', 'ObjVal_h2']].min(axis=1)
results.loc[results['init'] == 'with', 'Gap_math'] = 100*(results.loc[results['init']== 'with', 'ObjVal_h'] - results.loc[results['init']== 'with', 'ObjVal'])/results.loc[results['init']== 'with', 'ObjVal_h']
results.loc[results['init'] == 'with', 'Gap_build'] = 100*(results.loc[results['init']== 'with', 'ObjVal_h2'] - results.loc[results['init']== 'with', 'ObjVal_h'])/results.loc[results['init']== 'with', 'ObjVal_h2']
results.loc[results['init'] == 'with', 'Runtime_h'] = results.loc[results['init']== 'with', 'Runtime_h'].ffill(limit=100)
results.loc[results['init'] == 'with', 'Runtime_h2'] = results.loc[results['init']== 'with', 'Runtime_h2'].ffill(limit=100)

results.loc[results['init'] == 'with', 'Gap'] = results.loc[results['init']== 'with', 'Gap'].fillna(100)

results.to_csv('resultados/results_together.csv')

results_single = results.loc[results['approach'] == 'single']
results_multi = results.loc[results['approach'] == 'multi']

# pal_wL = sns.color_palette("Greys")

# pal_wL = ['#ededed', '#d1d1d1', '#adadad', '#828282', '#5c5c5c', '#2b2b2b']
# print(pal_wL.as_hex())


## Diferenciando por wL
tabla_comparadora = results.groupby(['n_N', 'wL', 'approach']).describe()[['Gap', 'Runtime']].round(2).reset_index()

results_wL = pd.DataFrame()
results_wL['approach'] = tabla_comparadora['approach']
results_wL['n_N'] = tabla_comparadora['n_N']
results_wL['wL'] = tabla_comparadora['wL']
results_wL['Gap'] = tabla_comparadora['Gap']['mean']
results_wL['Runtime'] = tabla_comparadora['Runtime']['mean']
results_wL['%Unsolved'] = np.round(100*(120 - tabla_comparadora['Gap']['count'])/120, 2)

results_wL = results_wL.pivot(index=['n_N','approach'], columns=['wL'])
# results_wL.to_excel('tex_files/tables/results_wL.xlsx')

## Diferenciando por initialisation
tabla_comparadora = results.groupby(['approach', 'n_N', 'init']).describe()[['Gap', 'Runtime', 'Gap_math', 'Gap_build', 'Runtime_h', 'Runtime_h2']].round(2).reset_index()

results_wL = pd.DataFrame()
results_wL['approach'] = tabla_comparadora['approach']
results_wL['n_N'] = tabla_comparadora['n_N']
results_wL['init'] = tabla_comparadora['init']
results_wL['Gap'] = tabla_comparadora['Gap']['mean']
results_wL['Gap_math'] = tabla_comparadora['Gap_math']['mean']
results_wL['Gap_build'] = tabla_comparadora['Gap_build']['mean']
results_wL['Runtime'] = tabla_comparadora['Runtime']['mean']
results_wL['Runtime_math'] = tabla_comparadora['Runtime_h']['mean']
results_wL['Runtime_build'] = tabla_comparadora['Runtime_h2']['mean']
results_wL['%Unsolved'] = np.round(100*(120 - tabla_comparadora['Gap']['count'])/120, 2)

results_wL = results_wL.pivot(index=['n_N','approach'], columns=['init'])
results_wL.to_excel('results_initialisation.xlsx')

## Diferenciando por porcentaje de barreras
tabla_comparadora = results.groupby(['approach', 'n_N', 'perc_B']).describe()[['Gap', 'Runtime']].round(2).reset_index()

results_wL = pd.DataFrame()
results_wL['approach'] = tabla_comparadora['approach']
results_wL['n_N'] = tabla_comparadora['n_N']
results_wL['perc_B'] = tabla_comparadora['perc_B']
results_wL['Gap'] = tabla_comparadora['Gap']['mean']
results_wL['Runtime'] = tabla_comparadora['Runtime']['mean']
results_wL['%Unsolved'] = np.round(100*(60 - tabla_comparadora['Gap']['count'])/60, 2)

results_wL = results_wL.pivot(index=['n_N','approach'], columns=['perc_B'])
# results_wL.to_excel('tex_files/tables/results_percB.xlsx')

# Diferenciando por k
tabla_comparadora = results.groupby(['approach', 'n_N', 'k']).describe()[['Gap', 'Runtime']].round(2).reset_index()

results_wL = pd.DataFrame()
results_wL['approach'] = tabla_comparadora['approach']
results_wL['n_N'] = tabla_comparadora['n_N']
results_wL['k'] = tabla_comparadora['k']
results_wL['Gap'] = tabla_comparadora['Gap']['mean']
results_wL['Runtime'] = tabla_comparadora['Runtime']['mean']
results_wL['%Unsolved'] = np.round(100*(80 - tabla_comparadora['Gap']['count'])/80, 2)

results_wL = results_wL.pivot(index=['n_N','approach'], columns=['k'])
# results_wL.to_excel('tex_files/tables/results_k.xlsx')

