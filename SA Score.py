import csv
import os
from BRSAScore import SAScorer
from rdkit import Chem
from rdkit.Chem import Draw

# Input and output paths
input_csv = "input_smiles.csv"
output_csv = "output_scores.csv"
contrib_dir = "contributions"

# Create scorer
scorer = SAScorer()

# Ensure output folder exists
os.makedirs(contrib_dir, exist_ok=True)

# Read input and calculate scores
results = []
with open(input_csv, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row["Name"]
        smiles = row["SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES for {name}: {smiles}")
            continue
        try:
            score, contrib = scorer.calculateScore(smiles)
            results.append({"Name": name, "SMILES": smiles, "BR-SAScore": round(score, 2)})
            
            # Save visualization
            img = Draw.MolToImage(mol, size=(300, 300))
            img.save(os.path.join(contrib_dir, f"{name}.png"))
        except Exception as e:
            print(f"Error scoring {name}: {e}")

# Write output CSV
with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Name", "SMILES", "BR-SAScore"])
    writer.writeheader()
    writer.writerows(results)

print(f"Done! Results saved to {output_csv} and images in {contrib_dir}")
