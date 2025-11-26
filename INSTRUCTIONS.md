# How to Run hgraph2graph (Antifungal Dataset)

Here are the steps to set up and run the training on your machine.

## 1. Clone the Repository

```bash
git clone https://github.com/jcatillo/hgraph2graph.git
cd hgraph2graph
```

## 2. Set up Environment

You need `conda` installed. Create a new environment and install dependencies:

```bash
# Create environment
conda create -n hgraph_env python=3.8 -y
conda activate hgraph_env

# Install RDKit (must be from conda-forge)
conda install -c conda-forge rdkit -y

# Install other dependencies
pip install torch networkx tqdm numpy
```

## 3. Preprocessing

We need to generate the vocabulary and preprocess the data into tensors.

```bash
# 1. Generate Vocabulary
# Note: We use 'all_single.txt' which has filtered out invalid multi-component molecules
python get_vocab.py --ncpu 8 < data/antifungal/all_single.txt > data/antifungal/vocab.txt

# 2. Preprocess Data
python preprocess.py --train data/antifungal/all_single.txt --vocab data/antifungal/vocab.txt --ncpu 8 --mode single

# 3. Organize Output
mkdir -p train_processed
mv tensor* train_processed/
```

## 4. Training

Now you can run the training. This code has been optimized to run on CPU if no GPU is detected.

```bash
# Create checkpoint directory
mkdir -p ckpt/antifungal-pretrained

# Run training
python train_generator.py --train train_processed/ --vocab data/antifungal/vocab.txt --save_dir ckpt/antifungal-pretrained
```

## Troubleshooting

- **"ModuleNotFoundError"**: Make sure you activated the environment (`conda activate hgraph_env`) and installed all requirements.
- **"ZeroDivisionError"**: This usually means the dataset is empty or too small. Ensure you are using `data/antifungal/all_single.txt`.
