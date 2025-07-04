# 🧪 Retrosynthesis and Molecule Scoring Tools

[![Python](https://img.shields.io/badge/python-3.7%2B-blue.svg)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-Installed-important)](https://www.rdkit.org/)
[![License: Research](https://img.shields.io/badge/license-research--only-lightgrey)]()

This repository provides a collection of Python tools for:

- Visualizing retrosynthesis trees from AiZynthFinder-style outputs  
- Scoring molecules for synthetic accessibility (BR-SAScore)  
- Evaluating synthetic complexity via SCScore models

---

## 📁 Contents

| Script           | Description                                                                 |
|------------------|-----------------------------------------------------------------------------|
| `json_to_word.py` | Converts retrosynthesis tree JSON into a structured Word document          |
| `SA Score.py`     | Calculates BR-SAScore from SMILES and saves both scores and structure images |
| `scscore.py`      | Computes synthetic complexity scores using SCScore deep learning models     |

---

## ⚙️ Requirements

Install required packages:

pip install numpy pandas rdkit python-docx anytree six
Additional needs:

BRSAScore.py module (for SA Score.py)

Trained SCScore model files (for scscore.py)

Proper input files: input_smiles.csv, output.json

📄 Usage
1. json_to_word.py
Input:

output.json: Retrosynthesis tree from AiZynthFinder

Output:

synthesis_tree.docx: Formatted tree with route scores

Run:

python json_to_word.py
2. SA Score.py
Input:

input_smiles.csv with Name and SMILES columns

Output:

output_scores.csv: Names, SMILES, BR-SAScores

contributions/: Molecule structure PNGs

Run:

python "SA Score.py"
Note: Requires BRSAScore.py in the same folder.

3. scscore.py
Input:

input_smiles.csv with SMILES column

SCScore models in models/ directory

Output:

SC_Score_output.csv: Scores from 3 models + average

Run:

python scscore.py
📂 Folder Structure
lua
.
├── json_to_word.py
├── SA Score.py
├── scscore.py
├── input_smiles.csv
├── output.json
├── models/
│   ├── full_reaxys_model_1024bool/
│   ├── full_reaxys_model_2048bool/
│   └── full_reaxys_model_1024uint8/
└── contributions/
✅ Example input_smiles.csv
csv
Name,SMILES
Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
Paracetamol,CC(=O)NC1=CC=C(C=C1)O
📌 Notes
json_to_word.py is useful for reports and visualizing retrosynthetic logic

SA Score.py combines BR-SAScore with images for synthesis planning

scscore.py applies fingerprint-based neural scoring for complexity

📜 License
This repository is intended for academic and research purposes only.
Please respect individual licenses of external tools such as RDKit, SCScore, and AiZynthFinder.


