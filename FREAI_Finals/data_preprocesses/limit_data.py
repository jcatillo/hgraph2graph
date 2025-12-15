import os
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Limit cleaned molecule data to specified number of SMILES')
parser.add_argument('--cleaned_path', type=str, required=True,
                    help='Path to the directory containing cleaned CSV files')
parser.add_argument('--output_path', type=str, required=True,
                    help='Path to the output directory for limited data')
parser.add_argument('--limit', type=int, default=680,
                    help='Number of SMILES to keep per file (default: 680)')
args = parser.parse_args()

# Paths to the cleaned data files
base_path = args.cleaned_path
output_path = args.output_path
files = ['cleaned_small.csv', 'cleaned_medium.csv', 'cleaned_large.csv']

# Create output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)
print(f"Output directory created: {output_path}\n")

# Limit to specified SMILES for each file (plus header)
limit = args.limit + 1

for file in files:
    file_path = os.path.join(base_path, file)
    output_file = os.path.join(output_path, file.replace('cleaned_', 'limited_'))
    print(f"Processing {file}...")
    
    # Read the CSV file
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    original_count = len(lines) - 1  # Subtract header
    print(f"  Original count: {original_count}")
    
    # Limit to specified lines (header + SMILES)
    limited_lines = lines[:limit]
    
    # Save to the limited directory
    with open(output_file, 'w') as f:
        f.writelines(limited_lines)
    
    print(f"  Limited to: {len(limited_lines) - 1} SMILES")
    print(f"  Saved to {output_file}\n")

print(f"All files have been limited to {args.limit} SMILES!")
