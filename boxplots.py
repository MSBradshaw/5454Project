import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

res_files = ['d_hi_res.csv','d_pkl_res.csv','d_reactome_res.csv','d_string_res.csv','rwr_hi_results.csv','rwr_pkl_results.csv','rwr_reactome_results.csv','rwr_string_results.csv']
df = None
for file in res_files:
    temp = pd.read_csv(file)
    temp['network'] = file.split('_')[1]
    temp['method'] = file.split('_')[0]
    if df is not None:
        df = pd.concat([temp,df])
    else:
        df = temp

# subset the columns
percent_df = df[['500 %','100 %', '50 %', '25 %', '10 %', 'disease', 'network','method']]

# make data tidy
tidy_percent_df = pd.melt(percent_df,['disease','network','method'], var_name='Top X Genes',value_name="value")
tidy_percent_df['Top X Genes'] = [int(x.replace(' %','') )for x in tidy_percent_df['Top X Genes']]
for m in set(tidy_percent_df['method']):
    plotting_df = tidy_percent_df[tidy_percent_df['method'] == m]
    splot = sns.boxplot(x="network", y="value", hue="Top X Genes", data=plotting_df, palette="Set1")
    plt.ylabel('Portion of genes recovered')
    plt.savefig('Figures/boxplot_'+m+'.png')
    plt.show()

percent_df_2 = percent_df[['500 %','100 %', '50 %', '25 %', '10 %','network']]

r = [1,2]
barWidth = 0.85

mean_series = percent_df_2.groupby('network').agg(pd.Series.mean)
p500 = np.array(list(mean_series['500 %']))
p100 = np.array(list(mean_series['100 %']))
p50 = np.array(list(mean_series['50 %']))
p25 = np.array(list(mean_series['25 %']))
p10 = np.array(list(mean_series['10 %']))

plt.bar(r, p500, color='#df7f20', edgecolor='white', width=barWidth, label='500')
plt.bar(r, p100, bottom=p500, color='#905998', edgecolor='white', width=barWidth, label='100')
plt.bar(r, p50, color='#59a257',bottom=p500+p100, edgecolor='white', width=barWidth, label='50')
plt.bar(r, p25, color='#477ca7',bottom=p500+p100+p50, edgecolor='white', width=barWidth, label='25')
plt.bar(r, p10, color='#cb3336',bottom=p500+p100+p50+p25, edgecolor='white', width=barWidth, label='10')
plt.xticks(r, ['PheKnowLater', 'String'])
plt.xlabel("Network")
plt.ylabel("Average portion of genes rediscovered")
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
plt.tight_layout()
plt.savefig('Figures/stackedbarplot.png')
plt.show()
