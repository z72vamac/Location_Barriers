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
fig, axes = plt.subplots(3, 2, sharex=True, figsize=(18, 8))
plt.subplots_adjust(wspace=0, hspace=0.2)

pal_wL = {0: "white", 50:"#828282"}
order_wL = [0, 50]

# Gap por wL
sns.boxplot(ax=axes[0, 0], x='n_N', y='Gap', hue='wL', hue_order= order_wL, data=results_single, legend=False, palette=pal_wL)
axes[0,0].set(title='single-commodity' ,xlabel = r"$|\mathcal{N}|$", ylabel = '\% Gap', ylim = (-1, 105))
axes[0,0].set_yticks([0, 20, 40, 60, 80, 100])

sns.boxplot(ax=axes[0, 1], x='n_N', y='Gap', hue='wL', hue_order= order_wL, data=results_multi, legend=False, palette=pal_wL)
axes[0,1].set(title='multi-commodity' ,xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (-1, 105))
axes[0,1].set_yticks([0, 20, 40, 60, 80, 100])
axes[0,1].tick_params(axis='y', labelleft=False)

# Runtime por wL
axes[1, 0].set(yscale='log')
sns.boxplot(ax=axes[1, 0], x='n_N', y='Runtime', hue='wL', hue_order= order_wL, data=results_single, legend=False, palette=pal_wL)
axes[1, 0].set(xlabel = r"$|\mathcal{N}|$", ylabel = 'Runtime (log scale)', ylim = (0.036, 6000))
axes[1, 0].set_yticks([0.036, 0.36, 3.6, 36, 360, 3600])
axes[1, 0].set_yticklabels([r"$3.6\times 10^{-2}$", r"$3.6\times 10^{-1}$", r"$3.6\times 10^{0}$", r"$3.6\times 10^{1}$", r"$3.6\times 10^{2}$", r"$3.6\times 10^{3}$"])


axes[1, 1].set(yscale='log')
sns.boxplot(ax=axes[1, 1], x='n_N', y='Runtime', hue='wL', hue_order= order_wL, data=results_multi, legend=False, palette=pal_wL)
axes[1, 1].set(xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (0.036, 6000))
axes[1, 1].set_yticks([0.036, 0.36, 3.6, 36, 360, 3600])
axes[1, 1].tick_params(axis='y', labelleft=False)

# Porcentaje que no encuentra solucion por wL
results['isnull'] = results.Gap.isnull()

prop_nan = results_single.groupby(['n_N', 'wL']).count().reset_index()
prop_nan['Gap'] = (120 - prop_nan['Gap'])/1.2
sns.barplot(ax=axes[2, 0], x='n_N', y='Gap', hue='wL', hue_order= order_wL, data=prop_nan, legend=False, linewidth=0.5, edgecolor='black', palette=pal_wL)
axes[2, 0].set(xlabel = r"$|\mathcal{N}|$", ylabel = '\% Problems without solution', ylim = (-1, 105))
axes[2, 0].set_yticks([0, 20, 40, 60, 80, 100])

prop_nan = results_multi.groupby(['n_N', 'wL']).count().reset_index()
prop_nan['Gap'] = (120 - prop_nan['Gap'])/1.2
sns.barplot(ax=axes[2, 1], x='n_N', y='Gap', hue='wL', hue_order= order_wL, data=prop_nan, legend=False, linewidth=0.5, edgecolor='black', palette=pal_wL)
axes[2, 1].set(xlabel = r"$|\mathcal{N}|$", ylim = (-1, 105), ylabel = None)
axes[2, 1].tick_params(axis='y', labelleft=False)
axes[2, 1].set_yticks([0, 20, 40, 60, 80, 100])


# Creando la leyenda por wL
white_patch = mpatches.Patch(facecolor='white', edgecolor='black', linewidth=0.5, label='0')
blue_patch = mpatches.Patch(facecolor='#828282', edgecolor='black', linewidth=0.5, label='50')

axes[0, 0].legend(handles=[white_patch, blue_patch], title = r'$\omega_L$', loc='upper left', frameon = False)

fig.align_ylabels(axes[:, 0])
# plt.tight_layout()

tikzplotlib.save('resultados/difference_between_wL.tex', encoding='utf-8')


## Diferenciando por initialisation
fig, axes = plt.subplots(3, 2, sharex=True, figsize=(18, 8))
plt.subplots_adjust(wspace=0, hspace=0.2)

pal_init = {'without': "white", 'with':"#828282"}
order_init = ['without', 'with']

# Gap por initialisation
sns.boxplot(ax=axes[0, 0], x='n_N', y='Gap', hue='init', hue_order= order_init, data=results_single, legend=False, palette=pal_init)
axes[0,0].set(title='single-commodity', xlabel = r"$|\mathcal{N}|$", ylabel = '\% Gap', ylim = (-1, 105))
axes[0,0].set_yticks([0, 20, 40, 60, 80, 100])

sns.boxplot(ax=axes[0, 1], x='n_N', y='Gap', hue='init', hue_order= order_init, data=results_multi, legend=False, palette=pal_init)
axes[0,1].set(title='multi-commodity', xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (-1, 105))
axes[0,1].set_yticks([0, 20, 40, 60, 80, 100])
axes[0,1].tick_params(axis='y', labelleft=False)

# Runtime por initialisation
axes[1, 0].set(yscale='log')
sns.boxplot(ax=axes[1, 0], x='n_N', y='Runtime', hue='init', hue_order= order_init, data=results_single, legend=False, palette=pal_init)
axes[1, 0].set(xlabel = r"$|\mathcal{N}|$", ylabel = 'Runtime', ylim = (0.036, 6000))
axes[1, 0].set_yticks([0.036, 0.36, 3.6, 36, 360, 3600])
axes[1, 0].set_yticklabels([r"$3.6\times 10^{-2}$", r"$3.6\times 10^{-1}$", r"$3.6\times 10^{0}$", r"$3.6\times 10^{1}$", r"$3.6\times 10^{2}$", r"$3.6\times 10^{3}$"])

axes[1, 1].set(yscale='log')
sns.boxplot(ax=axes[1, 1], x='n_N', y='Runtime', hue='init', hue_order= order_init, data=results_multi, legend=False, palette=pal_init)
axes[1, 1].set(xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (0.036, 6000))
axes[1, 1].set_yticks([0.036, 0.36, 3.6, 36, 360, 3600])
axes[1, 1].tick_params(axis='y', labelleft=False)

# Porcentaje que no encuentra solucion por initialisation
results['isnull'] = results.Gap.isnull()

prop_nan = results_single.groupby(['n_N', 'init']).count().reset_index()
prop_nan['Gap'] = (120 - prop_nan['Gap'])/1.2
sns.barplot(ax=axes[2, 0], x='n_N', y='Gap', hue='init', hue_order= order_init, data=prop_nan, legend=False, linewidth=0.5, edgecolor='black', palette=pal_init)
axes[2, 0].set(xlabel = r"$|\mathcal{N}|$", ylabel = '\% Problems without solution', ylim = (-1, 105))
axes[2, 0].set_yticks([0, 20, 40, 60, 80, 100])

prop_nan = results_multi.groupby(['n_N', 'init']).count().reset_index()
prop_nan['Gap'] = (120 - prop_nan['Gap'])/1.2
sns.barplot(ax=axes[2, 1], x='n_N', y='Gap', hue='init', hue_order= order_init, data=prop_nan, legend=False, linewidth=0.5, edgecolor='black', palette=pal_init)
axes[2, 1].set(xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (-1, 105))
axes[2, 1].tick_params(axis='y', labelleft=False)
axes[2, 1].set_yticks([0, 20, 40, 60, 80, 100])

# Creando la leyenda por initialisation
white_patch = mpatches.Patch(facecolor='white', edgecolor='black', linewidth=0.5, label='false')
blue_patch = mpatches.Patch(facecolor='#828282', edgecolor='black', linewidth=0.5, label='true')

axes[0, 0].legend(handles=[white_patch, blue_patch], title = 'initialisation', loc='upper left', frameon = False)

fig.align_ylabels(axes[:, 0])
# plt.tight_layout()

tikzplotlib.save('resultados/difference_between_initialisation.tex', encoding='utf-8')


## Diferenciando por porcentaje de barreras
fig, axes = plt.subplots(3, 2, sharex=True, figsize=(18, 8))
plt.subplots_adjust(wspace=0, hspace=0.2)

pal_perc = {'10.0%': "white", '20.0%': '#d1d1d1', '50.0%': '#adadad', '100%':"#828282"}
order_perc = ['10.0%', '20.0%', '50.0%', '100%']

## Gap por porcentaje de barreras
sns.boxplot(ax=axes[0, 0], x='n_N', y='Gap', hue='perc_B', hue_order=order_perc, data=results_single, legend=False, palette=pal_perc)
axes[0,0].set(title='single-commodity', xlabel = r"$|\mathcal{N}|$", ylabel = '\% Gap', ylim = (-1, 105))
axes[0,0].set_yticks([0, 20, 40, 60, 80, 100])

sns.boxplot(ax=axes[0, 1], x='n_N', y='Gap', hue='perc_B', hue_order=order_perc, data=results_multi, legend=False, palette=pal_perc)
axes[0,1].set(title='multi-commodity', xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (-1, 105))
axes[0,1].tick_params(axis='y', labelleft=False)

# Runtime por porcentaje de barreras
axes[1, 0].set(yscale='log')
sns.boxplot(ax=axes[1, 0], x='n_N', y='Runtime', hue='perc_B', hue_order=order_perc, data=results_single, legend=False, palette=pal_perc)
axes[1, 0].set_yticks([0.036, 0.36, 3.6, 36, 360, 3600])
axes[1, 0].set_yticklabels([r"$3.6\times 10^{-2}$", r"$3.6\times 10^{-1}$", r"$3.6\times 10^{0}$", r"$3.6\times 10^{1}$", r"$3.6\times 10^{2}$", r"$3.6\times 10^{3}$"])

axes[1, 1].set(yscale='log')
sns.boxplot(ax=axes[1, 1], x='n_N', y='Runtime', hue='perc_B', hue_order=order_perc, data=results_multi, legend=False, palette=pal_perc)
axes[1, 1].set(xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (0.036, 6000))
axes[1, 1].set_yticks([0.036, 0.36, 3.6, 36, 360, 3600])
axes[1, 1].tick_params(axis='y', labelleft=False)


# Porcentaje que no encuentra solucion por porcentaje de barreras
results['isnull'] = results.Gap.isnull()

prop_nan = results_single.groupby(['n_N', 'perc_B']).count().reset_index()
prop_nan['Gap'] = (60 - prop_nan['Gap'])/0.6
sns.barplot(ax=axes[2, 0], x='n_N', y='Gap', hue='perc_B', hue_order=order_perc, data=prop_nan, legend=False, linewidth=0.5, edgecolor='black', palette=pal_perc)
axes[2, 0].set(xlabel = r"$|\mathcal{N}|$", ylabel = '\% Problems without solution', ylim = (-1, 105))
axes[2, 0].set_yticks([0, 20, 40, 60, 80, 100])

prop_nan = results_multi.groupby(['n_N', 'perc_B']).count().reset_index()
prop_nan['Gap'] = (60 - prop_nan['Gap'])/0.6
sns.barplot(ax=axes[2, 1], x='n_N', y='Gap', hue='perc_B', hue_order=order_perc, data=prop_nan, legend=False, linewidth=0.5, edgecolor='black', palette=pal_perc)
axes[2, 1].set(xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (-1, 105))
axes[2, 1].tick_params(axis='y', labelleft=False)
axes[2, 1].set_yticks([0, 20, 40, 60, 80, 100])

# Creando la leyenda por initialisation
white_patch = mpatches.Patch(facecolor='white', edgecolor='black', linewidth=0.5, label=r'$10$')
dark1_patch = mpatches.Patch(facecolor='#d1d1d1', edgecolor='black', linewidth=0.5, label=r'$20$')
dark2_patch = mpatches.Patch(facecolor='#adadad', edgecolor='black', linewidth=0.5, label=r'$50$')
dark3_patch = mpatches.Patch(facecolor='#828282', edgecolor='black', linewidth=0.5, label=r'$100$')

axes[0, 0].legend(handles=[white_patch, dark1_patch, dark2_patch, dark3_patch], title = r'$\%|\mathcal B|$', loc='upper left', frameon = False)

fig.align_ylabels(axes[:, 0])
# plt.tight_layout()

tikzplotlib.save('resultados/difference_between_percentage.tex', encoding='utf-8')


# Diferenciando por k
fig, axes = plt.subplots(3, 2, sharex=True, figsize=(18, 8))
plt.subplots_adjust(wspace=0, hspace=0.2)

pal_k = {'1': "white", '10%': '#adadad', '25%':"#828282"}
order_k = ['1', '10%', '25%']

# Gap por k
sns.boxplot(ax=axes[0, 0], x='n_N', y='Gap', hue='k', hue_order=order_k, data=results_single, legend=False, palette=pal_k)
axes[0,0].set(title='single-commodity', xlabel = r"$|\mathcal{N}|$", ylabel = '\% Gap', ylim = (-1, 105))
axes[0,0].set_yticks([0, 20, 40, 60, 80, 100])

sns.boxplot(ax=axes[0, 1], x='n_N', y='Gap', hue='k', hue_order=order_k, data=results_multi, legend=False, palette=pal_k)
axes[0,1].set(title='multi-commodity', xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (-1, 105))
axes[0,1].tick_params(axis='y', labelleft=False)

# Runtime por k
axes[1, 0].set(yscale='log')
sns.boxplot(ax=axes[1, 0], x='n_N', y='Runtime', hue='k', hue_order=order_k, data=results_single, legend=False, palette=pal_k)
axes[1, 0].set(xlabel = r"$|\mathcal{N}|$", ylabel = 'Runtime', ylim = (0.036, 6000))
axes[1, 0].set_yticks([0.036, 0.36, 3.6, 36, 360, 3600])
axes[1, 0].set_yticklabels([r"$3.6\times 10^{-2}$", r"$3.6\times 10^{-1}$", r"$3.6\times 10^{0}$", r"$3.6\times 10^{1}$", r"$3.6\times 10^{2}$", r"$3.6\times 10^{3}$"])

axes[1, 1].set(yscale='log')
sns.boxplot(ax=axes[1, 1], x='n_N', y='Runtime', hue='k', hue_order=order_k, data=results_multi, legend=False, palette=pal_k)
axes[1, 1].set(xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (0.036, 6000))
axes[1, 1].set_yticks([0.036, 0.36, 3.6, 36, 360, 3600])
axes[1, 1].tick_params(axis='y', labelleft=False)


# Porcentaje que no encuentra solucion por k
results['isnull'] = results.Gap.isnull()
prop_nan = results_single.groupby(['n_N', 'k']).count().reset_index()
prop_nan['Gap'] = (80 - prop_nan['Gap'])/0.8
sns.barplot(ax=axes[2, 0], x='n_N', y='Gap', hue='k', hue_order=order_k, data=prop_nan, legend=False, linewidth=0.5, edgecolor='black', palette=pal_k)
axes[2, 0].set(xlabel = r"$|\mathcal{N}|$", ylabel = '\% Problems without solution', ylim = (-1, 105))
axes[2, 0].set_yticks([0, 20, 40, 60, 80, 100])

prop_nan = results_multi.groupby(['n_N', 'k']).count().reset_index()
prop_nan['Gap'] = (80 - prop_nan['Gap'])/0.8
sns.barplot(ax=axes[2, 1], x='n_N', y='Gap', hue='k', hue_order=order_k, data=prop_nan, legend=False, linewidth=0.5, edgecolor='black', palette=pal_k)
axes[2, 1].set(xlabel = r"$|\mathcal{N}|$", ylabel = None, ylim = (-1, 105))
axes[2, 1].tick_params(axis='y', labelleft=False)
axes[2, 1].set_yticks([0, 20, 40, 60, 80, 100])

# Creando la leyenda por initialisation
white_patch = mpatches.Patch(facecolor='white', edgecolor='black', linewidth=0.5, label=r'$1$')
dark2_patch = mpatches.Patch(facecolor='#adadad', edgecolor='black', linewidth=0.5, label=r'$10\%|\mathcal{N}|$')
dark3_patch = mpatches.Patch(facecolor='#828282', edgecolor='black', linewidth=0.5, label=r'$25\%|\mathcal{N}|$')

axes[0, 0].legend(handles=[white_patch, dark2_patch, dark3_patch], title = r'$k$', loc='upper left', frameon = False)

fig.align_ylabels(axes[:, 0])
# plt.tight_layout()

tikzplotlib.save('resultados/difference_between_k.tex', encoding='utf-8')


plt.show()
