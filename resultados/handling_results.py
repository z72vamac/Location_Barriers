import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


datos = pd.read_csv('resultados_buenos.csv')

datos = datos[datos['n_N'] <= 80]

tabla_comparadora = datos.groupby(['n_N', 'Lazy', 'A4']).describe()[['n_B', 'Gap', 'Runtime', 'ObjVal']].round(2).reset_index()

print(tabla_comparadora)

tabla_comparadora.to_excel('experiment.xlsx')