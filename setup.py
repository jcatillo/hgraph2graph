from setuptools import setup

setup(
    name="hgraph2graph",
    author="Wengong Jin",
    version="0.1.0",
    packages=['hgraph', 'props'],
    python_requires=">=3.6",
    install_requires=[
        "torch>=1.0.0",
        "networkx",
        "rdkit>=2019.03",
        "numpy",
        "pandas",
        "scikit-learn",
        "tqdm",
    ],
)
