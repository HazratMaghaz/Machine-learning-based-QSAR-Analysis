# Machine-learning-based-QSAR-Analysis
Machine learning based QSAR Model Generations, a berif guide step by step


# Atom-Based QSAR Model & Structural Optimization

## 🧪 Project Overview

This repository presents a complete **Atom-Based QSAR modeling workflow** focused on predicting IC50 values and suggesting structural modifications to improve compound potency.

> **Objective**: Build a predictive model for estimating IC50 values from molecular structures, then apply structural modifications (analogs) to propose activity-improving changes.

## 📦 Project Deliverables

- End-to-end code for descriptor calculation, model training, and SHAP interpretation.
- Model training notebooks (`RandomForest`, `XGBoost`, `ElasticNet`, `Ridge`).
- SHAP-based feature interpretation plots.
- Functional group optimization using analogs generated via RDKit.

## 🧬 Dataset (Client-owned, Not Included)

The dataset contains 16 small molecules with SMILES and experimental IC50 values (in µM). This data has **not been uploaded** due to client confidentiality.

## 🔬 Methodology Summary

### 1. **Descriptor Calculation**
- Used **Mordred** (Python) to compute >1600 molecular descriptors.
- Filtered to 1137 numeric descriptors after cleaning (constant, missing, non-numeric values removed).

### 2. **Modeling Approach**
- Trained multiple regression models:
  - `Random Forest`
  - `ElasticNet`
  - `XGBoost`
  - `Ridge Regression`

- Evaluated using:
  - **R² Score**
  - **MAE** (Mean Absolute Error)
  - **±0.5 pIC50 Accuracy**

### 3. **SHAP-Based Feature Analysis**
- Performed **SHAP beeswarm** and **waterfall** plots to identify top descriptors.
- Observed charge- and electrostatics-related descriptors (e.g., `PEOE_VSA7`, `MaxAbsPartialCharge`) as dominant.

### 4. **Analog Generation**
- Generated analogs of top 5 potent compounds using RDKit.
- Recalculated descriptors and evaluated with trained models.
- Compared IC50 shifts from parent to analog compounds.

| Compound | Parent pIC50 | Analog pIC50 | Δ Change |
|----------|--------------|--------------|----------|
| FI-3-8   | 4.92         | 4.82         | -0.10    |
| FI-3     | 4.88         | 4.75         | -0.13    |
| FI-3-1   | 4.77         | 4.69         | -0.08    |
| FI-3-12  | 4.79         | 4.70         | -0.09    |

---

## 📊 Final Model Performance

| Model       | R² Score | MAE    | Accuracy (±0.5) |
|-------------|----------|--------|-----------------|
| ElasticNet  | 0.7877   | 0.1004 | 100%            |
| RandomForest| Poor     | ~0.24  | 100% (overfit)  |
| XGBoost     | Poor     | ~0.29  | 75%             |

⚠️ Early models suffered from overfitting due to small sample size and high feature space. Final results improved after feature selection and analog design.

---

## 📈 Visualizations

| Plot Type         | Description                                |
|-------------------|--------------------------------------------|
| SHAP Beeswarm     | Shows global feature impact                |
| SHAP Waterfall    | Shows instance-level explanation           |
| IC50 Comparison   | Parent vs. analog activity                 |
| Feature Histogram | Distribution of key descriptors            |

---

## 🤖 Libraries Used

- `pandas`, `numpy`, `scikit-learn`
- `xgboost`, `matplotlib`, `seaborn`
- `mordred`, `rdkit`
- `shap`

---

## 💡 Conclusion

This project demonstrates how **explainable QSAR modeling** can guide rational compound optimization. Even with a limited dataset, using SHAP and analog design offers meaningful SAR (Structure-Activity Relationship) insights.

---

## 🔐 Confidentiality Notice

Certain folders (e.g., `DataSets/`, internal notes, raw CSVs) have been excluded for privacy. This project was developed as part of a client engagement and adheres to non-disclosure standards.

---

## 📫 Contact

If you have questions or wish to discuss collaboration or similar projects:

- **[Fiverr](https://www.fiverr.com/dromics05)**
- **[Email](maaz28608@gmail.com)**

---

