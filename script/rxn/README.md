# MCID Reaction Pipeline

Predicts metabolic reaction products from a compound database using SMIRKS-based reaction rules, then filters and annotates results with chemical properties and thermodynamic feasibility estimates.

## Input

| File | Format | Description |
|------|--------|-------------|
| `kegg_filtered.parquet` | Parquet | Substrate compounds. Required columns: `smiles`, `name`, `inchikey_main_block` |
| `rxn_info.csv` | CSV | Reaction rules. Required columns: `ReactionRule`, `Reaction`, `Description`, `ID` |

## Output

| File | Format | Description |
|------|--------|-------------|
| `kegg_1_rxn.parquet` | Parquet | Predicted products after filtering |

### Output columns

| Column | Description |
|--------|-------------|
| `substrate_smiles` | SMILES of the input substrate |
| `smiles` | SMILES of the predicted product |
| `reaction` | Reaction description |
| `reaction_rules` | SMIRKS pattern used |
| `reaction_id` | Reaction rule ID |
| `substrate_name` | Name from the input database |
| `inchi` | InChI of the product |
| `inchikey` | Full InChIKey of the product |
| `inchikey_first_block` | First block of InChIKey (connectivity layer) |
| `formula` | Molecular formula |
| `exact_mass` | Monoisotopic exact mass |
| `dg_estimate_kj_mol` | Estimated Gibbs free energy (kJ/mol) |
| `dg_uncertainty_kj_mol` | Uncertainty of the dG estimate |
| `thermodynamic_feasibility` | Feasibility category |
| `structural_changes` | Description of structural changes |

## Pipeline steps

1. Expand reaction rules (split multi-rule entries)
2. Apply reactions to each substrate (SMIRKS + sulfate conjugation + deoxyadenosine loss)
3. Remove invalid product SMILES
4. Remove products already present in the substrate database (by `inchikey_main_block`)
5. Estimate reaction thermodynamics (group contribution dG)
6. Remove thermodynamically unfavorable products

## Usage

```bash
cd script/rxn
python mcid_pipeline.py
```

## Dependencies

`rdkit`, `pandas`, `pyarrow`, `tqdm`, `matplotlib`, `Pillow`
