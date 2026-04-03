# sc-cell-state-benchmark

Benchmarking cell-state scoring methods in single-cell RNA-seq using biological perturbations, negative controls, and robustness tests to distinguish true biological signal from technical or stochastic noise.

## Overview

This project is a small, focused benchmarking study designed around a practical question in single-cell transcriptomics:

How reliable are common cell-state scoring methods, and how often do they capture real biology versus noise?

The project is motivated by a common challenge in single-cell RNA-seq analysis: automated cell-state classification and phenotype scoring can look convincing, but scores may be influenced by sequencing depth, batch effects, gene set composition, preprocessing choices, or stochastic variation. This repository evaluates how different scoring methods behave under known biological perturbation and under null or negative-control settings.

The main use case is immune activation scoring in single-cell data, with a focus on whether computational methods can robustly recover expected biology while resisting false signal.

## Project Goals

This repository aims to:

- Compare computational methods for quantifying phenotype activation across single-cell transcriptomes
- Establish negative-control and baseline frameworks for evaluating whether scores are biologically meaningful
- Test the sensitivity of scoring methods to technical perturbations such as normalization choices and downsampling
- Assess whether cell-state scores improve confidence in automated cell-state classification
- Provide interpretable summaries and visualizations suitable for scientific review

## Biological Rationale

In single-cell RNA-seq, each cell is represented by a high-dimensional gene expression profile. These profiles are used to infer both cell identity and cell state. However, expression measurements are noisy, sparse, and sensitive to technical effects. As a result, a high score for a pathway or phenotype does not automatically mean that the underlying biological program is truly active.

This project focuses on a setting where a biological signal is expected in advance: interferon response in peripheral blood immune cells. By choosing a dataset with a known perturbation, we can ask whether a scoring method detects the expected signal in the right cells and whether it remains stable under technical variation.

## Dataset Strategy

### Warm-Up Dataset

PBMC3k will be used as a small introductory dataset to practice:

- AnnData structure
- Preprocessing and filtering
- Normalization
- Clustering
- Marker gene interpretation
- Basic Scanpy workflow

This dataset is mainly for familiarization and pipeline validation.

### Main Benchmarking Dataset

The primary dataset for the project is the Kang PBMC interferon-beta stimulation dataset.

Why this dataset is suitable:

- It contains a known biological perturbation
- The expected signal is clear and biologically interpretable
- It allows direct comparison of stimulated versus unstimulated cells
- It provides a clean framework for benchmarking cell-state scoring approaches
- It is well suited for negative-control and robustness experiments

The main phenotype of interest will be interferon response or immune activation.

## Core Questions

The project is built around four main questions:

1. Do different cell-state scoring methods agree on which cells are activated?
2. Can these methods recover a known biological perturbation in the correct cell populations?
3. How often do they produce signal under negative-control conditions?
4. How stable are the scores when technical factors are changed?

## Specific Aims

### Aim 1: Compare Cell-State Scoring Methods

Implement and compare multiple approaches for quantifying phenotype activation in single cells.

Examples may include:

- Average expression of curated marker genes
- Scanpy module score
- Rank-based scoring
- Simple classifier-based score or probability

The first comparison will focus on interferon response gene sets.

### Aim 2: Build a Negative-Control Framework

Evaluate whether each scoring method produces inflated or misleading signal when no true biology should be present.

Negative controls may include:

- Random gene sets matched by size
- Shuffled gene identities
- Permuted condition labels
- Irrelevant pathway gene sets
- Depth-matched null comparisons

The goal is to determine whether a method is detecting meaningful biology or simply responding to generic data structure.

### Aim 3: Test Robustness to Technical Variation

Assess score stability under common analysis choices and perturbations.

Examples include:

- Different normalization strategies
- Different highly variable gene selections
- With and without batch correction
- Expression downsampling
- Varying gene set size

This aim evaluates whether conclusions are stable or fragile.

### Aim 4: Connect Scoring to Automated Cell-State Classification

Use the resulting scores to support or evaluate classification of cells into biological states such as stimulated versus unstimulated.

Questions include:

- Which scoring approach best separates known states
- Whether score-based classification is interpretable
- Whether a robust score improves confidence in cell-state assignment

## Planned Workflow

### 1. Data Loading and Preprocessing

- Load dataset into AnnData
- Perform quality control
- Normalize counts
- Identify highly variable genes
- Compute PCA, neighborhood graph, clustering, and UMAP
- Confirm expected cell populations and metadata structure

### 2. Biological Scoring

- Define interferon-response gene sets
- Compute multiple cell-state scores
- Compare distributions across conditions and cell types
- Visualize score patterns on UMAP and within annotated populations

### 3. Negative Controls

- Generate matched random gene sets
- Shuffle or permute labels where appropriate
- Compare real versus null score behavior
- Quantify false-positive tendencies

### 4. Robustness Experiments

- Downsample counts or cells
- Vary normalization and preprocessing choices
- Re-run scoring pipelines
- Summarize consistency of each method

### 5. Classification and Interpretation

- Classify stimulated versus unstimulated cells using scores
- Compare method performance across cell types
- Interpret which genes and pathways drive score behavior
- Summarize where methods are biologically plausible versus unreliable

## Expected Outputs

This project is intended to produce a compact but complete benchmarking analysis with:

- One reproducible analysis workflow
- Summary tables comparing scoring methods
- Negative-control benchmarks
- Robustness analyses
- Publication-style figures
- A concise technical summary of findings

## Example Figures

Planned visual outputs include:

- UMAP colored by cell type and condition
- Score distributions by cell type
- Real gene set scores versus matched random gene set scores
- Robustness plots under downsampling or preprocessing changes
- Classification performance comparison across scoring methods

## Methods and Tools

The project will primarily use Python-based single-cell analysis tools, including:

- Scanpy
- AnnData
- Pandas
- NumPy
- SciPy
- scikit-learn
- matplotlib

Potential extensions may include explainability or interpretability methods if classifier-based models are introduced later.

## Repository Structure

sc-cell-state-benchmark/
├── data/
├── notebooks/
├── src/
├── results/
├── figures/
├── README.md
└── requirements.txt

Suggested organization:

- data/: metadata, small processed files, or download instructions
- notebooks/: exploratory and analysis notebooks
- src/: reusable scoring and evaluation functions
- results/: summary tables and metrics
- figures/: final plots for reporting
- README.md: project overview and usage

## Initial Notebook Plan

1. 01_pbmc3k_scanpy_basics
   - Learn AnnData and standard preprocessing

2. 02_kang_data_preprocessing
   - Load and preprocess the perturbation dataset

3. 03_cell_state_scoring_methods
   - Implement and compare scoring approaches

4. 04_negative_controls
   - Evaluate null behavior and false positives

5. 05_robustness_analysis
   - Vary technical settings and measure stability

6. 06_classification_summary
   - Connect scores to automated state classification

## What This Project Demonstrates

This repository is designed to demonstrate the following skills:

- Understanding of single-cell transcriptomic analysis concepts
- Familiarity with Scanpy and AnnData workflows
- Benchmarking of computational methods
- Biological interpretation of pathway and phenotype activation
- Statistical thinking around baselines and negative controls
- Reproducible scientific programming and scientific communication

## Future Extensions

Possible future directions include:

- Validating on an additional public immune dataset
- Extending from interferon response to other immune phenotypes
- Adding explainable AI methods for classifier interpretation
- Comparing gene set scoring to latent embedding approaches
- Testing transferability across datasets and batches

## Status

Project in design and implementation phase.

## Author Motivation

This project was developed as a focused effort to strengthen hands-on experience in single-cell transcriptomics while building a method-oriented benchmark relevant to computational biology roles involving cell-state scoring, interpretability, and signal-versus-noise analysis.
