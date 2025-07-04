import pandas as pd
from rxnmapper import BatchedMapper

# --- File paths ---
input_csv = 'input_reactions.csv'       # ✅ Change to your actual input file
output_csv = 'mapped_output.csv'        # ✅ Output CSV file name

# --- Load data ---
df = pd.read_csv(input_csv, encoding='ISO-8859-1')

# --- Validate required columns ---
required_cols = {'R1', 'R2', 'P'}
if not required_cols.issubset(df.columns):
    raise ValueError(f"❌ Input file must contain columns: {required_cols}")

# --- Create reaction SMILES from R1, R2, P ---
def construct_reaction_smiles(row):
    r1 = row['R1'].strip() if pd.notna(row['R1']) else ''
    r2 = row['R2'].strip() if pd.notna(row['R2']) else ''
    p  = row['P'].strip()  if pd.notna(row['P'])  else ''
    
    if r1 and r2:
        return f"{r1}.{r2}>>{p}"
    elif r1:
        return f"{r1}>>{p}"
    else:
        return f">>{p}"  # Edge case: missing R1

df['Reaction_SMILES'] = df.apply(construct_reaction_smiles, axis=1)

# --- Initialize BatchedMapper ---
mapper = BatchedMapper(batch_size=32)

# --- Atom mapping (safe and robust) ---
results = list(mapper.map_reactions_with_info(df['Reaction_SMILES'].tolist()))

# --- Add mapping results ---
df['Mapped_Reaction'] = [r.get('mapped_rxn', '') if r else '' for r in results]
df['Mapping_Confidence'] = [r.get('confidence', 0.0) if r else 0.0 for r in results]

# --- Save to output CSV ---
df.to_csv(output_csv, index=False)

print(f"✅ Mapping completed successfully. Output saved to: {output_csv}")
