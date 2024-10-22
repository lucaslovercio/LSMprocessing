#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 15:28:28 2024

@author: lucas
"""

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC, LinearSVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import f1_score
from sklearn.utils import shuffle
import copy

plt.close('all')

# Load the CSV file
df = pd.read_csv('Proliferation_PC1-10_scores_UNICORNS.csv')

# Replace the 'genotype' column with binary labels: WT + Het as 0, Mut as 1
df['binary_genotype'] = df['genotype'].apply(lambda x: 1 if x == 'WT' else 0)

X = df[['PC1', 'PC2', 'PC3', 'PC4', 'PC5']] # WT vs Het-Mut

C_all = 1
# Binary target: Mut (1) vs WT + Het (0)
y = df['binary_genotype']

# Standardize the data
scaler = StandardScaler()
X_zscore = scaler.fit_transform(X)

# Convert X_zscore (numpy array) back to a DataFrame with the same columns and index as X
X_zscore_df = pd.DataFrame(X_zscore, columns=X.columns, index=X.index)

# Show the result
print(X_zscore_df.head())

X_zscore_df.to_csv('Proliferation_PC1-10_scores_UNICORNS_Zscore.csv')

# 1. Linear kernel with feature selection using L1 regularization
linear_svc = LinearSVC(penalty='l1', dual=False, max_iter=10000, C=C_all)  # 'dual=False' is required for L1 regularization
linear_svc.fit(X_zscore, y)

# Print coefficients for Linear SVM
print("Feature Coefficients for Linear SVM:")
for i, coef in enumerate(linear_svc.coef_[0]):
    print(f"PC{i + 1}: {coef:.4f}")

# Evaluate both models
print("Linear SVC (L1 regularization) Results:")
y_pred_linear = linear_svc.predict(X_zscore)
print(classification_report(y, y_pred_linear))

# Calculate the decision function (distance to the hyperplane) for each point in X_zscore
distances_to_hyperplane = linear_svc.decision_function(X_zscore)

# Convert y (Pandas Series) into a DataFrame
y_df = y.to_frame(name='genotype')
# Add the distances_to_hyperplane as a new column in y_df
y_df['distance_to_hyperplane'] = distances_to_hyperplane

print(y_df)

def permutation_test_f1(X, y, model, target_class=1, n_permutations=1000):
    # Fit the model to the original data and get predictions
    # model.fit(X, y)
    y_pred = model.predict(X)
    
    # Calculate the original F1-score for the Mut class
    original_f1 = f1_score(y, y_pred, pos_label=target_class)
    print('original_f1:' + str(original_f1))
    permuted_f1_scores = []

    for _ in range(n_permutations):
        y_permuted = shuffle(y)
        model.fit(X, y_permuted)
        y_permuted_pred = model.predict(X)
        
        # Calculate F1-score for the permuted labels
        permuted_f1 = f1_score(y_permuted, y_permuted_pred, pos_label=target_class)
        permuted_f1_scores.append(permuted_f1)

    # Calculate p-value
    p_value = np.sum(np.array(permuted_f1_scores) >= original_f1) / n_permutations
    return original_f1, p_value

# Assume 'original_model' is your trained model
copied_linear_svc = copy.deepcopy(linear_svc)
# Example usage
original_f1, p_value = permutation_test_f1(X_zscore, y, copied_linear_svc)
print(f'Original F1-Score for Mut class: {original_f1:.4f}')
print(f'P-Value from Permutation Test: {p_value:.4f}')


# ############# PLOTTING ################

# Create a scatter plot of PC4 vs PC5
plt.figure(figsize=(10, 6))

# Define colors for each genotype
colors = {'WT': 'blue', 'HET': 'green', 'MUT': 'red'}
# colors = {'WT': 'blue', 'MUT': 'red'}

# Plot each genotype with a different color
for genotype in colors:
    print(genotype)
    subset = df[df['genotype'] == genotype]
    # Subset the z-scored DataFrame using the index of the current genotype subset
    subset_zscore = X_zscore_df.loc[subset.index]
    plt.scatter(subset_zscore['PC3'], subset_zscore['PC4'], c=colors[genotype], label=genotype, s=100)

# Calculate the hyperplane
xlim = plt.xlim()
ylim = plt.ylim()
#xx = np.linspace(xlim[0], xlim[1])
x0 = xlim[0]
x1 = xlim[1]
y1 = (-linear_svc.coef_[0][2] * x1 - linear_svc.intercept_[0]) / linear_svc.coef_[0][3]  # WT vs Het-Mut
y0 = (-linear_svc.coef_[0][2] * x0 - linear_svc.intercept_[0]) / linear_svc.coef_[0][3]  # WT vs Het-Mut
# Plot the hyperplane as a dashed line
plt.plot([x0, x1], [y0, y1], 'k--', label='Hyperplane', linewidth=2)

# Customize the plot
plt.title('PC3 vs PC4 with Hyperplane of Separation')
plt.xlabel('PC3')
plt.ylabel('PC4')
plt.legend(title='Genotype')
plt.grid(True)
plt.xlim(xlim)
plt.ylim(ylim)

# Show the plot
plt.show()


from sklearn.model_selection import permutation_test_score

score_wt, perm_scores_wt, pvalue_wt = permutation_test_score(linear_svc, X_zscore, y, scoring='f1', n_permutations=1000)

print('pvalue_wt: ' + str(pvalue_wt))
