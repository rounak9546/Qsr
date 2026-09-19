# Scientific Methodology: QSAR & Computational Drug Discovery Workflow

This document details the complete 12-step scientific protocol implemented in this platform for Quantitative Structure-Activity Relationship (QSAR) modeling against Epidermal Growth Factor Receptor (EGFR).

---

## 1. Target Selection & Biological Rationale

- **Biological Target:** Epidermal Growth Factor Receptor (EGFR / ErbB-1).
- **Disease Indication:** Non-Small Cell Lung Cancer (NSCLC), glioblastoma multiforme, and colorectal cancer.
- **Mechanism of Action:** Competitive or covalent inhibition of the ATP-binding site located within the intracellular tyrosine kinase catalytic domain (residues 695–1022).
- **Oncogenic Rationale:** Activating mutations (exon 19 in-frame deletions, L858R point mutation) constitutively lock the kinase domain into the active state, driving uncontrolled cell survival through the PI3K-AKT and MAPK cascades.

---

## 2. Dataset Assembly & Preprocessing Protocol

1. **Ingestion:** Parse structured CSV records containing compound identifiers, chemical SMILES, and experimental bioassay measurements ($IC_{50}$).
2. **Duplicate Filtering:** Remove duplicate chemical structures identified by identical canonical SMILES strings or redundant ChEMBL registry identifiers.
3. **Invalid Record Removal:** Filter out compounds lacking valid chemical strings, unparsable SMILES, or missing/non-positive experimental activity records.
4. **Bioactivity Unit Normalization:** Convert reported potency measurements ($M$, $mM$, $\mu M$, $nM$, $pM$) into standard Molar concentration ($M$):
   $$\text{Molar Concentration } (M) = IC_{50} \times \text{Unit Factor}$$
5. **Logarithmic Transformation:** Convert molar inhibition concentrations into $pIC_{50}$:
   $$pIC_{50} = -\log_{10}(IC_{50} \text{ in Molar})$$

---

## 3. Chemoinformatics & Molecular Descriptors Calculation

Nine 1D/2D physicochemical and topological descriptors are computed directly from the molecular chemical graph:

1. **Molecular Weight (MW, Da):** Sum of standard atomic masses including implicit hydrogens.
2. **Calculated Octanol-Water Partition Coefficient ($\text{Log}P$):** Wildman-Crippen atom-additive contribution model.
3. **Hydrogen Bond Donors (HBD):** Count of nitrogen and oxygen atoms with attached hydrogens ($-OH, -NH_2, -NH-$).
4. **Hydrogen Bond Acceptors (HBA):** Count of all nitrogen and oxygen atoms with available lone pairs.
5. **Rotatable Bonds:** Count of non-ring, single bonds between non-terminal heavy atoms.
6. **Topological Polar Surface Area (TPSA, $\text{Å}^2$):** Sum of polar fragment contributions (Ertl et al., 2000).
7. **Ring Count:** Fundamental cycle basis of the molecular connection graph.
8. **Aromatic Ring Count:** Delocalized conjugated rings of 5 or 6 atoms.
9. **Fraction $sp^3$ ($\text{Fsp3}$):** Ratio of $sp^3$-hybridized carbons to total carbon count.

---

## 4. Machine Learning & QSAR Modeling Protocol

### Feature Standardization & Data Leakage Prevention
To prevent information leakage from test partitions into training routines:
1. Feature matrix $X$ is partitioned into training and test sets ($80/20$, $75/25$, or $70/30$) using a deterministic PRNG seed.
2. Mean ($\mu_j$) and standard deviation ($\sigma_j$) for each descriptor $j$ are calculated **strictly from the training partition**:
   $$\mu_j = \frac{1}{N_{\text{train}}} \sum_{i=1}^{N_{\text{train}}} X_{ij}, \quad \sigma_j = \sqrt{\frac{1}{N_{\text{train}}} \sum_{i=1}^{N_{\text{train}}} (X_{ij} - \mu_j)^2 + \epsilon}$$
3. Training and test matrices are standardized using these fixed training parameters:
   $$z_{ij} = \frac{X_{ij} - \mu_j}{\sigma_j}$$

### Algorithms Implemented
- **Multiple Linear Regression (OLS with Ridge Regularization):**
  $$\hat{\beta} = (X^T X + \lambda I)^{-1} X^T y$$
  where $\lambda = 10^{-4}$ ensures positive definiteness and prevents singular matrix errors.
- **Random Forest Regressor:**
  An ensemble of $B=20$ randomized decision trees trained using bootstrap aggregation (bagging). At each node split, a random subset of $m = \lfloor \sqrt{p} \rfloor$ features is evaluated to maximize variance reduction:
  $$\Delta \text{MSE} = \text{MSE}_{\text{parent}} - \left( \frac{N_L}{N}\text{MSE}_L + \frac{N_R}{N}\text{MSE}_R \right)$$
  Feature importance is computed as the total normalized variance reduction across all trees.

---

## 5. Model Validation & Statistical Evaluation

Models are evaluated on hold-out testing partitions and internal $K$-fold cross-validation:

- **Coefficient of Determination ($R^2$):**
  $$R^2 = 1 - \frac{\sum_{i=1}^n (y_i - \hat{y}_i)^2}{\sum_{i=1}^n (y_i - \bar{y})^2}$$
- **Root Mean Squared Error (RMSE):**
  $$\text{RMSE} = \sqrt{\frac{1}{n} \sum_{i=1}^n (y_i - \hat{y}_i)^2}$$
- **Mean Absolute Error (MAE):**
  $$\text{MAE} = \frac{1}{n} \sum_{i=1}^n |y_i - \hat{y}_i|$$
- **$K$-Fold Cross-Validation ($K=5, 10$):**
  $$\bar{R}^2_{CV} = \frac{1}{K} \sum_{k=1}^K R^2_k, \quad \text{SD}_{CV} = \sqrt{\frac{1}{K} \sum_{k=1}^K (R^2_k - \bar{R}^2_{CV})^2}$$

---

## 6. Applicability Domain (AD) Evaluation

Predictions are accompanied by an applicability domain check. For query molecule feature vector $x_{\text{query}}$:
1. Centroid $c$ of the standardized training dataset is derived.
2. The Euclidean distance $d_{\text{query}} = \|x_{\text{query}} - c\|_2$ is compared against the boundary threshold:
   $$\text{Threshold} = \bar{d}_{\text{train}} + 2.5 \cdot \text{SD}_{\text{train}}$$
3. If $d_{\text{query}} > \text{Threshold}$, the application alerts the researcher with an **extrapolation warning**.
