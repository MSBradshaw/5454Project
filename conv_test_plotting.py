import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from os import listdir

df = None

for file in listdir('./'):
    if 'converging_test' in file:
        print(file)
        temp = pd.read_csv(file)
        if df is None:
            df = temp
        else:
            df = pd.concat([df,temp])
df = df.reset_index()
totals = []
for i in range(df.shape[0]):
    print(i)
    totals.append(df[['500 count','100 count','50 count','25 count','10 count']].iloc[i,:].sum())
df['totals'] = totals

sns.scatterplot(data=df,x='threshold',y='totals')
plt.xlim([0.000009,0.000015])
plt.savefig('threshold_v_total_rediscovered.png')
plt.show()

df['hours'] = df['time'] / (60 * 60)
sns.scatterplot(data=df,x='threshold',y='hours')
plt.xlim([0.000009,0.000015])
plt.savefig('threshold_v_time.png')
plt.show()




df = None

for file in listdir('./Weighted_Converging_Results'):
    if 'converging_test' in file:
        print(file)
        temp = pd.read_csv('./Weighted_Converging_Results/' + file)
        if df is None:
            df = temp
        else:
            df = pd.concat([df,temp])
df = df.reset_index()
totals = []
for i in range(df.shape[0]):
    print(i)
    totals.append(df[['500 count','100 count','50 count','25 count','10 count']].iloc[i,:].sum())

df['totals'] = totals

sns.scatterplot(data=df,x='threshold',y='totals')
plt.xlim([0.000009,0.000015])
plt.savefig('weighted_threshold_v_total_rediscovered.png')
plt.show()

df['hours'] = df['time'] / (60 * 60)
sns.scatterplot(data=df,x='threshold',y='hours')
plt.xlim([0.000009,0.000015])
plt.savefig('weighted_threshold_v_time.png')
plt.show()
