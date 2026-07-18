# -*- coding: utf-8 -*-

@author: HP
"""


import pandas as pd
import numpy as np
import lightgbm as lgb
import shap
from sklearn.model_selection import RepeatedKFold
from sklearn.metrics import r2_score
from scipy import stats

# 数据预处理 ----
dff = pd.read_stata('D:/真菌分析/data_pre/2015_mmsepredict_conti.dta') 

dff.set_index('SampleID', inplace=True)
data = dff.replace('.', np.nan)

# 提取特征和目标变量
target = data['MMSE_score'].astype('float64')
feature = data.drop(columns='MMSE_score')

X = feature.values
y = target.values

# 设置迭代参数 ----
n_repeats = 100 
n_splits = 10    # 10 折交叉验证
rkf = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=2026)

# 用于存储结果的容器
all_r2_scores = []
all_pearson_r = []
# 用于存储特征重要性 (SHAP)
all_shap_values = np.zeros(X.shape) 

print(f"开始 {n_repeats} 次重复的 {n_splits} 折交叉验证...")

# 3. 核心循环 ----
# RepeatedKFold 会生成 n_repeats * n_splits 次迭代
for i, (train_index, test_index) in enumerate(rkf.split(X)):
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    
    # 训练模型
    model = lgb.LGBMRegressor(
        objective='regression',
        learning_rate=0.05,
        n_estimators=100,  #
        importance_type='gain',
        verbosity=-1
    )
    model.fit(X_train, y_train)
    
    # 预测并计算当前 fold 的得分
    y_pred = model.predict(X_test)
    
    # 记录 R2
    r2 = r2_score(y_test, y_pred)
    all_r2_scores.append(r2)
    
    # 记录 Pearson r
    r, _ = stats.pearsonr(y_test, y_pred)
    all_pearson_r.append(r)
    explainer = shap.TreeExplainer(model)
    # 仅针对测试集计算 SHAP 值，累加到全局
    all_shap_values[test_index] += explainer.shap_values(X_test) / n_repeats

    if (i + 1) % 50 == 0:
        print(f"进度: 已完成 {i + 1} / {n_repeats * n_splits} 个 Fold")

# 5. 结果整理与稳健统计 ----

# 将结果转换为每轮重复的平均值（每 10 个 fold 为一轮）
r2_per_repeat = [np.mean(all_r2_scores[i:i+10]) for i in range(0, len(all_r2_scores), 10)]
pearson_per_repeat = [np.mean(all_pearson_r[i:i+10]) for i in range(0, len(all_pearson_r), 10)]

print("\n" + "="*40)
print("── 最终稳健性统计 (基于 {} 次重复) ──".format(n_repeats))
print(f"R2 Score:      中位数 = {np.median(r2_per_repeat):.4f}, [IQR: {np.percentile(r2_per_repeat, 25):.4f} - {np.percentile(r2_per_repeat, 75):.4f}]")
print(f"Pearson r:     中位数 = {np.median(pearson_per_repeat):.4f}, [IQR: {np.percentile(pearson_per_repeat, 25):.4f} - {np.percentile(pearson_per_repeat, 75):.4f}]")
print(f"R2 标准误(SEM): {np.std(r2_per_repeat) / np.sqrt(n_repeats):.4f}")
print("="*40)

# 特征重要性汇总 ----
importance_df = pd.DataFrame({
    'feature': feature.columns,
    'mean_abs_shap': np.abs(all_shap_values).mean(axis=0)
}).sort_values(by='mean_abs_shap', ascending=False)

print(importance_df.head(10))


import os

# 设置保存路径
save_path = r'D:\真菌aging\投稿\NA revise\machine learning'

# 检查文件夹是否存在，如果不存在则创建
if not os.path.exists(save_path):
    os.makedirs(save_path)
    print(f"创建文件夹: {save_path}")

# 导出特征重要性 (importance_df)
# 包含特征名称和平均绝对 SHAP 值
importance_file = os.path.join(save_path, 'feature_importance_shap.csv')
importance_df.to_csv(importance_file, index=False, encoding='utf-8-sig')

# 导出每轮重复的 Pearson r 结果 (pearson_per_repeat)
# 将列表转换为 DataFrame 以便导出
pearson_df = pd.DataFrame({
    'repeat_round': range(1, n_repeats + 1),
    'mean_pearson_r': pearson_per_repeat,
    'mean_r2_score': r2_per_repeat  # 顺带把 R2 也存了，方便后续分析
})
pearson_file = os.path.join(save_path, 'model_performance_repeats.csv')
pearson_df.to_csv(pearson_file, index=False, encoding='utf-8-sig')

print("\n" + "-"*30)
print(f"数据已成功导出至:\n1. {importance_file}\n2. {pearson_file}")
print("-"*30)


#-----------------------shuffle the label to validate----------

# 数据准备 (保持原样)
dff = pd.read_stata('D:/真菌分析/data_pre/2015_mmsepredict_conti.dta')
dff.set_index('SampleID', inplace=True)
data = dff.replace('.', np.nan)
target = data['MMSE_score'].astype('float64')
feature = data.drop(columns='MMSE_score')
X = feature.values
y = target.values

# 设置参数：平衡严谨性与速度
n_splits = 10
n_permutations = 100 
kf = KFold(n_splits=n_splits, shuffle=True, random_state=2026) # 

def run_fast_evaluation(X_input, y_input):
    """单次 10 折交叉验证"""
    r_list = []
    # 直接在 X 上循环，不使用 Repeat
    for train_index, test_index in kf.split(X_input):
        X_train, X_test = X_input[train_index], X_input[test_index]
        y_train, y_test = y_input[train_index], y_input[test_index]
        
        model = lgb.LGBMRegressor(n_estimators=100, learning_rate=0.05, verbosity=-1)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        
        r, _ = stats.pearsonr(y_test, y_pred)
        r_list.append(r)
    return np.mean(r_list)

# 计算真实观察值
print("正在计算真实 Pearson r...")
obs_r = run_fast_evaluation(X, y)
print(f"Observed Pearson r: {obs_r:.4f}")

#执行置换检验 (100 次打乱)
print(f"开始执行 {n_permutations} 次快速置换检验...")
null_r_distribution = []

for p in range(n_permutations):
    # 随机打乱 y
    y_shuffled = np.random.permutation(y)
    
    # 每次打乱跑一次完整的 10 折 CV
    p_r = run_fast_evaluation(X, y_shuffled)
    null_r_distribution.append(p_r)
    
    if (p + 1) % 10 == 0:
        print(f"进度: {p + 1}% 已完成")

# 结果汇总
p_permutation = (np.sum(np.array(null_r_distribution) >= obs_r) + 1) / (n_permutations + 1)

# 可视化
plt.figure(figsize=(8, 5))
sns.histplot(null_r_distribution, kde=True, color="skyblue", label='Null Distribution')
plt.axvline(obs_r, color='red', linestyle='--', label=f'Observed r (p < {p_permutation:.3f})')
plt.title('Permutation Test (n=100)')
plt.xlabel('Pearson correlation (r)')
plt.legend()
plt.tight_layout()
plt.show()

print(f"分析完成！最终 P 值范围: p < {p_permutation:.3f}")
