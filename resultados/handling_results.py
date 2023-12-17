import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
# from scipy import stats
import matplotlib
import tikzplotlib
import os
import matplotlib.patches as mpatches

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

results.rename(columns={'single':'approach'}, inplace=True)

results['approach'] = results['approach'].astype('category')
results['approach'] = results['approach'].replace({True:'single', False:'multi'})
results['approach'].cat.reorder_categories(['single', 'multi'])

# results['index'] = np.arange(len(results))

results['Gap'] *= 100
# results['ObjVal_h2'] = results[['ObjVal_h', 'ObjVal_h2']].min(axis=1)
# results.loc[results['init'] == True, 'Runtime'] = results.loc[results['init']== True, 'Runtime'].fillna(3600)
vector = ['Gap' , 'Runtime','NodeCount','ObjVal','Runtime_h','Runtime_h2','ObjVal_h','ObjVal_h2']

results.loc[results['init'] == True, vector] = results.loc[results['init'] == True].groupby(['n_N', 'perc_B', 'k', 'approach', 'wL','lazy','A4']).transform(lambda x: x.fillna(x.mean())).reset_index().loc[:, vector]

# results.loc[results['init'] == True, 'NodeCount'] = results.loc[results['init']== True, 'NodeCount'].fillna(method='ffill', limit=100)

# results.loc[results['init'] == True, 'ObjVal_h'] = results.loc[results['init']== True, 'ObjVal_h'].ffill(limit=100)
# results.loc[results['init'] == True, 'Runtime_h'] = results.loc[results['init']== True, 'Runtime_h'].fillna(method='ffill', limit=100)
# results.loc[results['init'] == True, 'Runtime_h2'] = results.loc[results['init']== True, 'Runtime_h2'].fillna(method='ffill', limit=100)

# results.loc[results['init'] == True, 'Gap'] = results.loc[results['init']== True, 'Gap'].fillna(100)

results.to_csv('resultados/results_together.csv')

sns.set(style="darkgrid")

matplotlib.rcParams["font.sans-serif"] = 'Times New Roman'
matplotlib.rcParams["font.serif"] = 'Times New Roman'
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['text.usetex'] = True

g = sns.catplot(x='n_N', y='Gap', hue='approach', col='wL', kind='box',  hue_order=['single', 'multi'], legend=False, data=results, sharey=True, palette="Blues")

# sns.move_legend(g, "lower right", bbox_to_anchor=(-2, -1))

tikzplotlib.save('resultados/difference_between_wL(gap).tex', encoding='utf-8')

g.set(xlabel = r"$|\mathcal{N}|$",
       ylabel = '\% Gap',
       ylim = (0, 100))

# ax.set_yticklabels(np.arange(0, 100, 10))
# g._legend.set_title("Approach")

for ax, title in zip(g.axes.flat, [r"$\omega_{L}=0$", r"$\omega_{L}=50$"]):
    ax.set_title(title)

lightblue_patch = mpatches.Patch(color='#B2CDDE', label='Single-commodity')
blue_patch = mpatches.Patch(color='#4884AF', label='Multi-commodity')

plt.legend(handles=[lightblue_patch, blue_patch], title = 'approach', loc='upper right', frameon = False)

plt.tight_layout()
plt.show()
