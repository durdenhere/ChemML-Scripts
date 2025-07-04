import math, sys, random, os
import numpy as np
import time
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import json
import gzip
import six
import pandas as pd

# Project path
project_root = os.path.dirname(os.path.abspath(__file__))

# Define constants
score_scale = 5.0
min_separation = 0.25
FP_len = 1024
FP_rad = 2

def sigmoid(x):
    return 1 / (1 + math.exp(-x))

class SCScorer():
    def __init__(self, score_scale=score_scale):
        self.vars = []
        self.score_scale = score_scale
        self._restored = False

    def restore(self, weight_path, FP_rad=FP_rad, FP_len=FP_len):
        self.FP_len = FP_len
        self.FP_rad = FP_rad
        self._load_vars(weight_path)
        print('Restored variables from {}'.format(weight_path))

        if 'uint8' in weight_path or 'counts' in weight_path:
            def mol_to_fp(self, mol):
                if mol is None:
                    return np.zeros((self.FP_len,), dtype=np.uint8)
                fp = AllChem.GetMorganFingerprint(mol, self.FP_rad, useChirality=True)
                fp_folded = np.zeros((self.FP_len,), dtype=np.uint8)
                for k, v in six.iteritems(fp.GetNonzeroElements()):
                    fp_folded[k % self.FP_len] += v
                return np.array(fp_folded)
        else:
            def mol_to_fp(self, mol):
                if mol is None:
                    return np.zeros((self.FP_len,), dtype=np.float32)
                return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, self.FP_rad, nBits=self.FP_len,
                    useChirality=True), dtype=np.bool)
        self.mol_to_fp = mol_to_fp

        self._restored = True
        return self

    def smi_to_fp(self, smi):
        if not smi:
            return np.zeros((self.FP_len,), dtype=np.float32)
        return self.mol_to_fp(self, Chem.MolFromSmiles(smi))

    def apply(self, x):
        if not self._restored:
            raise ValueError('Must restore model weights!')
        # Each pair of vars is a weight and bias term
        for i in range(0, len(self.vars), 2):
            last_layer = (i == len(self.vars)-2)
            W = self.vars[i]
            b = self.vars[i+1]
            x = np.matmul(x, W) + b
            if not last_layer:
                x = x * (x > 0)  # ReLU
        x = 1 + (score_scale - 1) * sigmoid(x)
        return x

    def get_score_from_smi(self, smi='', v=False):
        if not smi:
            return ('', 0.)
        
        # Check if SMILES is valid
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            if v: print(f'Invalid SMILES: {smi}')
            return ('', 0.)
        
        fp = np.array((self.smi_to_fp(smi)), dtype=np.float32)
        if sum(fp) == 0:
            if v: print('Could not get fingerprint')
            cur_score = 0.
        else:
            # Run
            cur_score = self.apply(fp)
            if v: print('Score: {}'.format(cur_score))
        
        # Get canonical SMILES
        smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True)
        return (smi, cur_score)

    def _load_vars(self, weight_path):
        if weight_path.endswith('pickle'):
            import pickle
            with open(weight_path, 'rb') as fid:
                self.vars = pickle.load(fid)
                self.vars = [x.tolist() for x in self.vars]
        elif weight_path.endswith('json.gz'):
            with gzip.GzipFile(weight_path, 'r') as fin:
                json_bytes = fin.read()
                json_str = json_bytes.decode('utf-8')
                self.vars = json.loads(json_str)
                self.vars = [np.array(x) for x in self.vars]


def load_all_models(models_dir):
    """Load all three SCScore models"""
    models = {}
    
    # Model 1: 1024bool
    print("Loading model 1024bool...")
    model_1024bool = SCScorer()
    model_1024bool.restore(
        os.path.join(models_dir, 'full_reaxys_model_1024bool', 'model.ckpt-10654.as_numpy.json.gz'),
        FP_len=1024
    )
    models['1024bool'] = model_1024bool
    
    # Model 2: 2048bool
    print("Loading model 2048bool...")
    model_2048bool = SCScorer()
    model_2048bool.restore(
        os.path.join(models_dir, 'full_reaxys_model_2048bool', 'model.ckpt-10654.as_numpy.json.gz'),
        FP_len=2048
    )
    models['2048bool'] = model_2048bool
    
    # Model 3: 1024uint8
    print("Loading model 1024uint8...")
    model_1024uint8 = SCScorer()
    model_1024uint8.restore(
        os.path.join(models_dir, 'full_reaxys_model_1024uint8', 'model.ckpt-10654.as_numpy.json.gz'),
        FP_len=1024
    )
    models['1024uint8'] = model_1024uint8
    
    return models


def score_smiles_with_all_models(smi, models, verbose=False):
    """Score a single SMILES with all three models"""
    results = {
        'SMILES': '',
        'Score_1024bool': 0.0,
        'Score_2048bool': 0.0,
        'Score_1024uint8': 0.0,
        'Avg_SCScore': 0.0
    }
    
    # Check if SMILES is valid first
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        if verbose:
            print(f"Skipping invalid SMILES: {smi}")
        return None
    
    scores = []
    canonical_smi = ''
    
    for model_name, model in models.items():
        try:
            canon_smi, score = model.get_score_from_smi(smi, v=verbose)
            if canon_smi:  # Only proceed if we got a valid canonical SMILES
                canonical_smi = canon_smi
                results[f'Score_{model_name}'] = score
                scores.append(score)
                if verbose:
                    print(f"{model_name}: {score:.4f}")
            else:
                if verbose:
                    print(f"Could not process SMILES with {model_name}: {smi}")
                return None
        except Exception as e:
            if verbose:
                print(f"Error processing {smi} with {model_name}: {str(e)}")
            return None
    
    if len(scores) == 3:  # All models succeeded
        results['SMILES'] = canonical_smi
        results['Avg_SCScore'] = np.mean(scores)
        return results
    else:
        return None


if __name__ == '__main__':
    input_csv = 'input_smiles.csv'
    output_csv = 'SC_Score_output.csv'
    
    # Load input SMILES
    try:
        df = pd.read_csv(input_csv)
        if 'SMILES' not in df.columns:
            raise ValueError("Input CSV must have a 'SMILES' column")
        print(f"Loaded {len(df)} SMILES from {input_csv}")
    except Exception as e:
        print(f"Error loading input CSV: {str(e)}")
        sys.exit(1)
    
    # Load all models
    models_dir = os.path.join(project_root, 'models')
    try:
        models = load_all_models(models_dir)
        print("All models loaded successfully!")
    except Exception as e:
        print(f"Error loading models: {str(e)}")
        sys.exit(1)
    
    # Process all SMILES
    results = []
    valid_count = 0
    invalid_count = 0
    
    print(f"\nProcessing {len(df)} SMILES...")
    for idx, smi in enumerate(df['SMILES']):
        if idx % 100 == 0:
            print(f"Processed {idx}/{len(df)} SMILES...")
        
        result = score_smiles_with_all_models(smi, models, verbose=False)
        if result is not None:
            results.append(result)
            valid_count += 1
        else:
            invalid_count += 1
    
    # Write output
    if results:
        out_df = pd.DataFrame(results)
        out_df.to_csv(output_csv, index=False)
        print(f"\nResults saved to {output_csv}")
        print(f"Valid SMILES processed: {valid_count}")
        print(f"Invalid SMILES skipped: {invalid_count}")
        print(f"Success rate: {valid_count/(valid_count+invalid_count)*100:.1f}%")
        
        # Show sample results
        print("\nSample results:")
        print(out_df.head())
    else:
        print("No valid SMILES were processed!")
