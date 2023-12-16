import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
# from scipy import stats
import matplotlib
import tikzplotlib
import os

# cwd = os.getcwd()
# files = os.listdir(cwd) 
# print("Files in %r: %s" % (cwd, files))

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

# results['index'] = np.arange(len(results))

results['Gap'] *= 100
# results['ObjVal_h2'] = results[['ObjVal_h', 'ObjVal_h2']].min(axis=1)
# results.loc[results['init'] == True, 'Runtime'] = results.loc[results['init']== True, 'Runtime'].fillna(3600)
vector = ['Gap' , 'Runtime','NodeCount','ObjVal','Runtime_h','Runtime_h2','ObjVal_h','ObjVal_h2']

results.loc[results['init'] == True, vector] = results.loc[results['init'] == True].groupby(['n_N', 'perc_B', 'k', 'single', 'wL','lazy','A4']).transform(lambda x: x.fillna(x.mean())).reset_index().loc[:, vector]

# results.loc[results['init'] == True, 'NodeCount'] = results.loc[results['init']== True, 'NodeCount'].fillna(method='ffill', limit=100)

# results.loc[results['init'] == True, 'ObjVal_h'] = results.loc[results['init']== True, 'ObjVal_h'].ffill(limit=100)
# results.loc[results['init'] == True, 'Runtime_h'] = results.loc[results['init']== True, 'Runtime_h'].fillna(method='ffill', limit=100)
# results.loc[results['init'] == True, 'Runtime_h2'] = results.loc[results['init']== True, 'Runtime_h2'].fillna(method='ffill', limit=100)

# results.loc[results['init'] == True, 'Gap'] = results.loc[results['init']== True, 'Gap'].fillna(100)

results.to_csv('resultados/results_together.csv')

g = sns.catplot(x='n_N', y='Gap', hue='single', col='wL', kind='box', data=results, aspect=1, sharey=True)

sns.set(style="darkgrid")


matplotlib.rcParams['axes.unicode_minus'] = False

tikzplotlib.save('tex_files/figures/difference_between_wL.tex', encoding='utf-8')

plt.show()
