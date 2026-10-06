import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("benchmark_results.csv")

df_constants = df[df['config'].str.startswith('WT_') & (df['N'] == 10**6)].copy()

plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
sns.lineplot(data=df_constants, x='config', y='build_time_ms', hue='distribution', marker='o')
plt.title('Время Insert (Build Time) vs B_leaf')
plt.xticks(rotation=45)
plt.ylabel('Время (мс)')
plt.grid()

plt.subplot(1, 2, 2)
sns.lineplot(data=df_constants, x='config', y='query_time_ms', hue='distribution', marker='o')
plt.title('Время Aggregate (Query Time) vs B_leaf')
plt.xticks(rotation=45)
plt.ylabel('Время (мс)')
plt.grid()
plt.tight_layout()
plt.show()

df_baselines = df[df['config'].isin(['WT_Balanced(256_16)', 'Sqrt_Decomposition', 'Pruned_Treap'])].copy()

plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
sns.lineplot(data=df_baselines[df_baselines['distribution'] == 'uniform'],
             x='N', y='query_time_ms', hue='config', marker='o')
plt.xscale('log'); plt.yscale('log')
plt.title('Время запросов vs N (Uniform)')
plt.ylabel('Время (мс)'); plt.grid()

plt.subplot(1, 2, 2)
sns.lineplot(data=df_baselines[df_baselines['distribution'] == 'zipf'],
             x='N', y='query_time_ms', hue='config', marker='o')
plt.xscale('log'); plt.yscale('log')
plt.title('Время запросов vs N (Zipf)')
plt.ylabel('Время (мс)'); plt.grid()
plt.tight_layout()
plt.show()

df_mem = df[df['N'] == 10**6].copy()
df_mem['memory_mb'] = df_mem['size_in_bytes'] / (1024 * 1024)

plt.figure(figsize=(10, 6))
sns.barplot(data=df_mem, x='config', y='memory_mb', hue='distribution')
plt.title('Потребление памяти [N=10^6]')
plt.xticks(rotation=45)
plt.ylabel('Память (МБ)')
plt.grid(axis='y')
plt.show()