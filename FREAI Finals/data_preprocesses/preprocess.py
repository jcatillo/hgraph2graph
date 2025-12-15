# main.py
import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem


def remove_fragmented_smiles(smiles_list):
    cleaned = []
    for smi in smiles_list:
        if smi is None:
            continue
        smi = str(smi).strip()
        if not smi:
            continue
        if "." in smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            cleaned.append(smi)
    return cleaned


def canonicalize_smiles(smiles_list):
    canonicalize_smiles = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            can_smi = Chem.MolToSmiles(mol, canonical=True)
            canonicalize_smiles.append(can_smi)
    return canonicalize_smiles

def deduplicate_smiles(smiles_list, keep_order=True):
    if keep_order:
        # preserves order
        return list(dict.fromkeys(smiles_list))
    return list(set(smiles_list))


import argparse
from pathlib import Path
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Clean SMILES from a CSV column.")
    p.add_argument("-i", "--input", required=True, help="Input CSV file path.")
    p.add_argument("-o", "--output", required=True, help="Output file path (csv).")
    p.add_argument("--smiles-col", default="SMILES", help="Name of the SMILES column.")
    p.add_argument(
        "--write-smiles-only",
        action="store_true",
        help="Output CSV will contain only the (processed) SMILES column.",
    )

    p.add_argument("--remove-fragments", action="store_true", help="Drop SMILES containing '.'")
    p.add_argument("--canonicalize", action="store_true", help="Canonicalize SMILES with RDKit")
    p.add_argument("--dedup", action="store_true", help="Deduplicate SMILES")

    p.add_argument(
        "--no-keep-order",
        action="store_true",
        help="If set, dedup may reorder SMILES (uses set).",
    )
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input)

    if args.smiles_col not in df.columns:
        raise ValueError(
            f"SMILES column '{args.smiles_col}' not found. Columns: {list(df.columns)}"
        )

    smiles = df[args.smiles_col].tolist()
    input_n = len(smiles)

    # ✅ run only what user asked for, in a fixed order
    if args.remove_fragments:
        smiles = remove_fragmented_smiles(smiles)

    if args.canonicalize:
        smiles = canonicalize_smiles(smiles)

    if args.dedup:
        smiles = deduplicate_smiles(smiles, keep_order=not args.no_keep_order)

    # output
    if args.write_smiles_only:
        out_df = pd.DataFrame({args.smiles_col: smiles})
    else:
        # If you removed/deduped, row count may change; SMILES-only is safest.
        out_df = pd.DataFrame({args.smiles_col: smiles})

    out_path = Path(args.output)
    out_df.to_csv(out_path, index=False)

    print("Done.")
    print(f"Input rows:  {input_n}")
    print(f"Output rows: {len(smiles)}")
    print(f"Wrote: {out_path}")

if __name__ == "__main__":
    main()

