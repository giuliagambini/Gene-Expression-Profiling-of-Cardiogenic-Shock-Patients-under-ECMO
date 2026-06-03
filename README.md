# Gene Expression Profiling of Cardiogenic Shock Patients under ECMO

[![R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📌 Project Overview
This repository contains the code, analysis report, and presentation for the multi-omics analysis of cardiogenic shock patients supported by Extracorporeal Membrane Oxygenation (ECMO). The primary goal is to identify a robust transcriptomic predictive signature for early clinical outcomes and map the underlying systemic biological drivers.

## 📂 Repository Structure
* **`/code`**: R scripts used for data pre-processing, machine learning classification, and functional integration.
* **`/report`**: The full academic report detailing the methodology, results, and biological interpretation (PDF format).
* **`/presentation`**: The slide deck used for the final project discussion.
* **`/graphs`**: Key visualizations generated during the analysis (e.g., LASSO shrinkage, STRING networks, rCCA plots).

## 🧰 Methodology & Tools
The analytical pipeline integrates machine learning with network-based biology:
1. **Pre-processing & Exploratory Analysis**: PCA, K-means, Hierarchical Clustering.
2. **Supervised Machine Learning (Feature Selection & Classification)**: 
   * Random Forest (RF)
   * Linear Discriminant Analysis (LDA)
   * LASSO Regularization
   * SCUDO (Signature-based Clustering for Diagnostic Outcomes)
3. **Multi-omics & Functional Integration**:
   * **g:Profiler** & **STRING**: Gene Ontology enrichment and Protein-Protein Interaction networks.
   * **pathfindR**: Active subnetwork identification and pathway clustering.
   * **rCCA (regularized Canonical Correlation Analysis)**: Integration of transcriptomic profiles with clinical parameters.

## 📊 Key Findings
* **Predictive Power**: Coordinated pathway analysis drastically outperforms isolated markers in predicting ECMO outcomes.
* **Biological Drivers**: The outcome is heavily driven by the systemic dysregulation of lipid homeostasis and inflammatory networks.
* **Clinical Relevance**: Network-based functional integration provides a solid, data-driven framework for patient stratification in intensive care settings.

## 🚀 How to Run the Code
1. Clone this repository:
```bash
   git clone [https://github.com/giuliagambini/](https://github.com/giuliagambini/)[Gene-Expression-Profiling-of-Cardiogenic-Shock-Patients-under-ECMO].git

## 👤 Author
**Giulia Gambini** *University of Trento* 📧 Email: giulia.gambini@studenti.unitn.it  
🔗 LinkedIn: [Giulia Gambini](https://www.linkedin.com/in/giulia-gambini-778644299)
