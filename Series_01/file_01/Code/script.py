# =======================================
#  QSAR Pipeline (Parent + Analog) Code
#  From SMILES to Predictions
# =======================================

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions
from mordred import Calculator, descriptors
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from sklearn.pipeline import make_pipeline

# ================================
# Step 1: Load Initial Dataset
# ================================
df = pd.read_csv("FI_Series_QSAR_Data_with_pIC50.csv")  # Contains Compound_ID, SMILES, IC50_µM, pIC50

# Clean and check molecules
df["mol"] = df["SMILES"].apply(Chem.MolFromSmiles)
df = df[df["mol"].notnull()].reset_index(drop=True)

# ================================
# Step 2: Generate Analogs using RDKit Reaction Templates
# ================================
reactions = [
    "[#6:1][CH3:2]>>[#6:1][OH]",      # Replace CH3 with OH
    "[#6:1][H]>>[#6:1][Cl]",          # Add Cl to carbon
    "[#6:1][H]>>[#6:1][F]"            # Add F to carbon
]

analogs = []
for idx, row in df.iterrows():
    mol = row['mol']
    for smirks in reactions:
        try:
            rxn = rdChemReactions.ReactionFromSmarts(smirks)
            products = rxn.RunReactants((mol,))
            for product in products:
                smi = Chem.MolToSmiles(product[0])
                analogs.append({
                    "Parent_ID": row["Compound_ID"],
                    "Parent_SMILES": row["SMILES"],
                    "Analog_SMILES": smi
                })
        except:
            continue

analog_df = pd.DataFrame(analogs)
analog_df.drop_duplicates(subset=["Analog_SMILES"], inplace=True)
analog_df["mol"] = analog_df["Analog_SMILES"].apply(Chem.MolFromSmiles)
analog_df = analog_df[analog_df["mol"].notnull()].reset_index(drop=True)
analog_df["Type"] = "Analog"

# ================================
# Step 3: Calculate Mordred Descriptors
# ================================
calc = Calculator(descriptors, ignore_3D=True)

# Parent compounds
desc_parent = calc.pandas(df["mol"])
parent_df = pd.concat([df[["Compound_ID", "SMILES", "pIC50"]].reset_index(drop=True), desc_parent], axis=1)
parent_df["Type"] = "Parent"

# Analogs
desc_analogs = calc.pandas(analog_df["mol"])
analog_ready = pd.concat([analog_df[["Parent_ID", "Parent_SMILES", "Analog_SMILES"]].reset_index(drop=True), desc_analogs], axis=1)

# ================================
# Step 4: Merge and Clean Final Dataset
# ================================
combined = pd.concat([parent_df, analog_ready], ignore_index=True)
X = combined.select_dtypes(include=[np.number]).drop(columns=["pIC50"], errors="ignore")
y = combined["pIC50"] if "pIC50" in combined.columns else None

X = X.replace([np.inf, -np.inf], np.nan).dropna(axis=1)
X = X.loc[:, X.nunique() > 1]

# Train-test split only on parent compounds
X_train = X[combined["Type"] == "Parent"]
y_train = y[combined["Type"] == "Parent"]
X_analog = X[combined["Type"] == "Analog"]

# ================================
# Step 5: Model Training & Prediction
# ================================
model = make_pipeline(StandardScaler(), ElasticNet(alpha=0.1, l1_ratio=0.5, random_state=42))
model.fit(X_train, y_train)

# Predict pIC50 for analogs
y_pred = model.predict(X_analog)
combined.loc[combined["Type"] == "Analog", "Predicted_pIC50"] = y_pred

# Save Final Output
combined.to_csv("QSAR_Parent_Analog_Predictions.csv", index=False)
print("✅ Done: Saved prediction results to QSAR_Parent_Analog_Predictions.csv")
