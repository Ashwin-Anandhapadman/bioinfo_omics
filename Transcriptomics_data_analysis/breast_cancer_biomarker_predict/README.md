# Breast cancer prediction using Machine learning 


This folder contains the code written by me to perform transcriptomics data analysis on clinical RNA seq data obtained from TCGA database in this research paper:https://doi.org/10.1186/s12859-022-04965-8.

I have to tried to replicate the results shown in the above paper using my data analysis workflow. This ML pipeline finds the best set of features from the dataset and identifies potential transcriptomic biomarkers for breast cancer.

## More info on the dataset and approach:
In this study, the dataset was obtained from 762 breast cancer patients and 138 normal subjects with solid tissue. The approach involved the application of various machine learning (ML) algorithms, organized into three main groups:

### 1. Feature Selection: 
Compared four feature selection methods to identify the most significant features:
1.1.ANOVA
1.2.Mutual Information
1.3. Extra Trees Classifier
1.4. Logistic Regression (LGR)

### 2. Feature Extraction: 
Utilized Principal Component Analysis (PCA) as our feature extraction technique.

### 3. Classification Algorithms: 
A total of 13 classification algorithms were implemented:
3.1. Logistic Regression (LGR)
3.2. Support Vector Machine (SVM)
3.3. Bagging
3.4. Gaussian Naive Bayes
3.5. Decision Tree
3.6. Gradient Boosting Decision Tree
3.7. K Nearest Neighbors (KNN)
3.8. Bernoulli Naive Bayes
3.9. Random Forest
3.10. AdaBoost
3.11. Extra Trees
3.12. Linear Discriminant Analysis (LDA)
3.13. Multilayer Perceptron (MLP) 
