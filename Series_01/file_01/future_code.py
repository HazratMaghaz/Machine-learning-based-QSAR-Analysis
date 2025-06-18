# 📁 Step 1: Load Dataset and Prepare SMILES
import pandas as pd
from rdkit import Chem
from mordred import Calculator, descriptors
from tqdm import tqdm

# Load original dataset
original_df = pd.read_csv("FI_Series_QSAR_Data_with_pIC50.csv")
original_df = original_df[['Compound_ID', 'SMILES', 'pIC50']].dropna()
original_df['mol'] = original_df['SMILES'].apply(Chem.MolFromSmiles)

# Filter valid molecules
original_df = original_df[original_df['mol'].notnull()].reset_index(drop=True)

# Add molecule type
original_df['Type'] = 'Parent'

# View shape
print(f"✅ Loaded original dataset: {original_df.shape}")

# 📌 We'll use this DataFrame to build the model later.


# 📁 Step 2: Generate Analogs
from rdkit.Chem import AllChem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

# Select top compounds (highest pIC50)
top_parents = original_df.sort_values(by='pIC50', ascending=False).head(3)

# Define a SMARTS-based reaction for -Br to -CH3 substitution
reaction = ReactionFromSmarts('[c:1][Br:2]>>[c:1][CH3:2]')

# Store generated analogs
generated = []

for idx, row in top_parents.iterrows():
    mol = row['mol']
    products = reaction.RunReactants((mol,))
    for p in products:
        try:
            smi = Chem.MolToSmiles(p[0])
            generated.append({
                'Parent_ID': row['Compound_ID'],
                'Parent_SMILES': row['SMILES'],
                'Analog_SMILES': smi,
                'Type': 'Analog'
            })
        except:
            continue

analogs_df = pd.DataFrame(generated).drop_duplicates('Analog_SMILES')
print(f"✅ Total analogs generated: {analogs_df.shape[0]}")


# 📁 Step 3: Calculate Mordred Descriptors for All
# Merge both datasets
analog_mols = [Chem.MolFromSmiles(smi) for smi in analogs_df['Analog_SMILES']]
analogs_df['mol'] = analog_mols
analogs_df = analogs_df[analogs_df['mol'].notnull()].reset_index(drop=True)

# Combine both
all_mols_df = pd.concat([
    original_df[['Compound_ID', 'SMILES', 'pIC50', 'mol', 'Type']],
    analogs_df.rename(columns={'Parent_ID': 'Compound_ID', 'Analog_SMILES': 'SMILES'})
], ignore_index=True)

# Setup Mordred
calc = Calculator(descriptors, ignore_3D=True)

# Calculate descriptors
print("🧪 Calculating Mordred descriptors...")
desc_df = calc.pandas(tqdm(all_mols_df['mol']))

# Merge back
final_df = pd.concat([all_mols_df[['Compound_ID', 'SMILES', 'pIC50', 'Type']], desc_df], axis=1)

# Clean up
final_df = final_df.select_dtypes(include=['number'])
final_df.dropna(axis=1, inplace=True)
print(f"✅ Final descriptor shape: {final_df.shape}")


# 📁 Step 4: Train Model on Parents Only
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np

# Split parent data only
parent_df = final_df[final_df['Type'] == 'Parent']
X = parent_df.drop(columns=['pIC50'])
X = X.loc[:, X.nunique() > 1]
y = parent_df['pIC50']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train ElasticNet model
model = make_pipeline(StandardScaler(), ElasticNet(alpha=0.1, l1_ratio=0.5, random_state=42))
model.fit(X_train, y_train)

# Evaluate
y_pred = model.predict(X_test)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
r2 = r2_score(y_test, y_pred)
print(f"✅ RMSE: {rmse:.4f}")
print(f"✅ R² Score: {r2:.4f}")


# 📁 Step 5: Predict Activity for Analogs
analog_df = final_df[final_df['Type'] == 'Analog']
X_analogs = analog_df.drop(columns=['pIC50'])
X_analogs = X_analogs[X_train.columns]

# Predict pIC50 for analogs
analog_df['Predicted_pIC50'] = model.predict(X_analogs)
print("✅ Prediction complete for analogs")

# Save everything
final_df.to_csv("QSAR_Mordred_Parent_Analog_Full.csv", index=False)
analog_df.to_csv("QSAR_Analog_Predictions.csv", index=False)
print("✅ Files saved: All descriptors and predictions")