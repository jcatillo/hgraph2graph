from sklearn.model_selection import train_test_split
import os
import random
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Split limited molecule data into train/test sets and combine them')
parser.add_argument('--limited_path', type=str, required=True, 
                    help='Path to the directory containing limited CSV files')
parser.add_argument('--output_path', type=str, required=True,
                    help='Path to the output directory for split data')
parser.add_argument('--test_size', type=float, default=0.2,
                    help='Proportion of test set (default: 0.2)')
parser.add_argument('--seed', type=int, default=42,
                    help='Random seed for reproducibility (default: 42)')
args = parser.parse_args()

# Paths
limited_path = args.limited_path
output_path = args.output_path
files = ['limited_small.csv', 'limited_medium.csv', 'limited_large.csv']

# Set random seed
random.seed(args.seed)

# Create output directory
os.makedirs(output_path, exist_ok=True)
print(f"Output directory created: {output_path}\n")

# Split each file into 80:20 train-test
for file in files:
    file_path = os.path.join(limited_path, file)
    print(f"Processing {file}...")
    
    # Read the CSV file
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Separate header and data
    header = lines[0]
    data = lines[1:]
    
    print(f"  Total SMILES: {len(data)}")
    
    # Split 80:20
    train_data, test_data = train_test_split(data, test_size=args.test_size, random_state=args.seed)
    
    print(f"  Train: {len(train_data)}, Test: {len(test_data)}")
    
    # Create output filenames
    base_name = file.replace('limited_', '').replace('.csv', '')
    train_file = os.path.join(output_path, f'{base_name}_train.csv')
    test_file = os.path.join(output_path, f'{base_name}_test.csv')
    
    # Save train set
    with open(train_file, 'w') as f:
        f.write(header)
        f.writelines(train_data)
    
    # Save test set
    with open(test_file, 'w') as f:
        f.write(header)
        f.writelines(test_data)
    
    print(f"  Saved: {train_file}")
    print(f"  Saved: {test_file}\n")

print("All files have been split into 80:20 train-test sets!")

# Combine all training data
print("\nCombining all training data...")
all_train_smiles = []
for file in files:
    base_name = file.replace('limited_', '').replace('.csv', '')
    train_file = os.path.join(output_path, f'{base_name}_train.csv')
    with open(train_file, 'r') as f:
        lines = f.readlines()[1:]  # Skip header
        all_train_smiles.extend(lines)

# Shuffle the combined training data
random.shuffle(all_train_smiles)

# Save combined training data
combined_train_file = os.path.join(output_path, 'train.txt')
with open(combined_train_file, 'w') as f:
    f.writelines(all_train_smiles)
print(f"Combined training data: {len(all_train_smiles)} SMILES")
print(f"Saved to: {combined_train_file}")

# Combine all test data
print("\nCombining all test data...")
all_test_smiles = []
for file in files:
    base_name = file.replace('limited_', '').replace('.csv', '')
    test_file = os.path.join(output_path, f'{base_name}_test.csv')
    with open(test_file, 'r') as f:
        lines = f.readlines()[1:]  # Skip header
        all_test_smiles.extend(lines)

# Shuffle the combined test data
random.shuffle(all_test_smiles)

# Save combined test data
combined_test_file = os.path.join(output_path, 'test.txt')
with open(combined_test_file, 'w') as f:
    f.writelines(all_test_smiles)
print(f"Combined test data: {len(all_test_smiles)} SMILES")
print(f"Saved to: {combined_test_file}")

print("\nAll data has been combined and shuffled!")

