## Hemoglobin Complex Analysis Project

### Overview
This project analyzes the structural properties and interactions of the hemoglobin complex within the human proteome. The focus is on identifying key components that stabilize the hemoglobin tetramer, such as salt bridges and hotspots, and understanding the interactions between protein monomers and heme groups. Additionally, the connectivity of hemoglobin subunits within the human proteome is explored to gain insights into their roles in various biological pathways.

### Repository Contents
- **`scripts/`**: Contains Python scripts used for the analysis.
- **`data/`**: Includes input data files such as PDB and DSSP files.
- **`README.md`**: This readme file.

### Getting Started

#### Prerequisites
- Python 3.7.2 or later
- Required Python libraries: Numpy, NetworkX
- UCSF Chimera for 3D visualization
- Unix-based system for running commands

#### Installation
1. Clone the repository:
   ```sh
   git clone https://github.com/yourusername/hemoglobin-complex-analysis.git
   cd hemoglobin-complex-analysis
   ```
2. Install required Python libraries:
   ```sh
   pip install numpy networkx
   ```

### Usage

#### Data Preparation
1. Download the hemoglobin PDB structure (ID: 1GZX) from the PDB website and place it in the `data/` directory.
2. Generate DSSP files for the full hemoglobin, all possible trimers, and each monomer using the DSSP program. Place these files in the `data/` directory.

#### Running the Analysis
1. **PDB Parsing**: Run the script to parse the PDB file and analyse interactions between residues:
   ```sh
   python scripts/parsepdb.py    #HEME or MONO
   ```
2. **DSSP Parsing**: Run the script to analyse the accessible surface area and identify interacting residues:
   ```sh
   python scripts/parsedssp.py   #RASA or ASA
   ```
3. **Protein-Protein Interaction Network**: Run the script to analyse the protein-protein interaction network:
   ```sh
   python scripts/network.py
   ```

### Visualization
Use UCSF Chimera to visualize the 3D structures and interactions. Load the PDB files and view the interactions highlighted in the analysis.


