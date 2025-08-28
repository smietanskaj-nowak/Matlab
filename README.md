# Matlab Scripts for Protein Structure Analysis and Visualization

This repository contains a collection of **MATLAB scripts** designed for the analysis and visualization of protein structures.  
The scripts are intended for researchers, students, and bioinformaticians who work with structural biology data, such as PDB (Protein Data Bank) files or molecular dynamics output.

---

## Features

- Import and process protein structure data
- Compute structural parameters (bond lengths, bond angles, torsions)
- Visualize protein backbones and side-chain connectivity
- Generate customizable 2D/3D plots of protein geometry
- Tools for comparative analysis across multiple protein conformations

---

## Contents

- **Data files** – example structural data for testing (e.g., connectivity matrices, bond parameters)  
- **Analysis scripts** – MATLAB functions for:
  - Loading protein structure data - **load_structure_data.m**
  - Calculating distances, bond angles, and torsional angles - **tables_coefficients**
  - Extracting geometric descriptors - **calc_basevectors**
- **Visualization scripts** – tools for:
  - Rendering protein backbones - **plot_structure.m**
  - Highlighting residues, bonds, or regions of interest - **select_mol**
  - Generating plots for publication or reports - **cylinder2P**

---

## Getting Started

### Requirements
- MATLAB R2020a or newer (earlier versions may work but are untested)
- Bioinformatics Toolbox (recommended for working with PDB files)

### Usage
1. Clone or download this repository:
   ```bash
   git clone https://github.com/your-username/matlab-protein-structure-tools.git
