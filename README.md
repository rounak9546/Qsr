# Integration of QSAR, Bioinformatics and Machine Learning for Cancer Drug Discovery (EGFR)

An interactive, client-side computational biology and chemoinformatics platform demonstrating the end-to-end drug discovery pipeline: from molecular structure and descriptors to Quantitative Structure-Activity Relationship (QSAR) machine learning modeling, statistical validation, and candidate prioritization for **EGFR (Epidermal Growth Factor Receptor)**.

Developed for academic presentation, B.Tech Biotechnology / Bioinformatics degree defense, portfolio demonstration, and open-source computational drug discovery research.

---

## 🌟 Platform Highlights

- **100% Client-Side & Static Host Ready**: Zero backend or localhost requirement. Deployable directly to Vercel, GitHub Pages, or Netlify.
- **Authentic Target & Benchmark Data**: Curated dataset of 120 authentic EGFR small-molecule inhibitors from ChEMBL (Target: `CHEMBL203`), including approved drugs (Gefitinib, Erlotinib, Osimertinib, Afatinib, Lapatinib).
- **Chemoinformatics Engine**: Computes 9 essential 1D/2D molecular descriptors in-browser (MW, LogP, HBD, HBA, Rotatable Bonds, TPSA, Ring Count, Aromatic Rings, Fraction Csp3) from SMILES with Lipinski Rule of 5 and Veber rule evaluation.
- **Pure JavaScript Machine Learning**:
  - Multiple Linear Regression (Ordinary Least Squares with L2 regularization)
  - Ridge Regression ($\lambda = 0.01$)
  - Random Forest Regressor (Decision Tree ensemble with bootstrap bagging and feature sub-sampling)
  - K-Fold Cross-Validation ($k=5, 10$)
  - Authentic metrics: $R^2$, RMSE, MAE, MSE
  - Data leakage prevention (feature scaling fitted strictly on training partition)
- **Interactive Visualizations**:
  - Actual vs. Predicted scatter plot with $y = x$ identity line
  - Residual error distribution plot
  - Molecular descriptor histograms
  - 2D Chemical space PCA (Principal Component Analysis)
  - Pearson correlation matrix heatmap
  - Feature importance ranking bar chart
  - Model benchmark comparison chart
- **Candidate Prioritization & Virtual Screening**: High-throughput ranking of candidate molecules by predicted $pIC_{50}$ with CSV export.
- **Interactive 2D Structure Drawer**: Visualizes molecular graphs on HTML5 Canvas.

---

## 📂 Project Structure

```
/qsr
│
├── index.html                  # Main SPA containing all 12 navigation sections
├── vercel.json                 # Static routing & caching headers for Vercel
├── .gitignore                  # Git ignore rules
│
├── css/
│   └── style.css               # Dark scientific dashboard theme, responsive styling
│
├── js/
│   ├── app.js                  # Master application orchestrator & lifecycle manager
│   ├── ui.js                   # Navigation, toast notifications, modals, spinners, tables
│   ├── data.js                 # Dataset loading, CSV parser, preprocessor, summary stats, PCA
│   ├── descriptors.js          # Molecular descriptor calculation & 2D chemical structure canvas
│   ├── qsar.js                 # QSAR transformations (IC50 to pIC50, unit conversions, AD domain)
│   ├── ml.js                   # Pure JS ML models (Linear, Ridge, Random Forest, CV, metrics)
│   ├── charts.js               # Scientific visualizations using Chart.js with canvas fallback
│   └── validation.js           # CSV schema verification, SMILES validation, numeric range checks
│
├── data/
│   └── egfr_activity.csv       # Curated authentic dataset of verified EGFR inhibitors from ChEMBL
│
├── assets/
│   ├── images/
│   │   └── egfr_pathway.svg    # EGFR signaling pathway and kinase inhibition diagram
│   └── icons/
│       └── favicon.svg         # Scientific chemical hexagon icon
│
├── README.md                   # Complete documentation, setup, deployment guide
├── DATA_SOURCES.md             # Dataset provenance, ChEMBL / PubChem references
├── METHODOLOGY.md              # Detailed mathematical formulas & scientific protocol
└── VIVA.md                     # 25+ academic defense questions & answers
```

---

## 🚀 Quick Start & Local Setup

Because the platform is purely static, you can run it locally with any simple HTTP server:

### Option 1: Python HTTP Server
```bash
# In the project root directory:
python -m http.server 8000
```
Open your browser at `http://localhost:8000`.

### Option 2: Node.js (npx serve)
```bash
npx serve .
```

### Option 3: VS Code Live Server
Right-click `index.html` and select **"Open with Live Server"**.

---

## 🌐 Deploying to Vercel

The repository includes a preconfigured `vercel.json`. You can deploy in seconds:

1. Push your repository to GitHub:
   ```bash
   git init
   git add .
   git commit -m "Initial commit of EGFR QSAR Bioinformatics Platform"
   git remote add origin https://github.com/your-username/egfr-qsar-discovery.git
   git push -u origin main
   ```
2. Go to [Vercel](https://vercel.com/) and click **"Add New Project"**.
3. Import your GitHub repository.
4. Framework Preset: select **"Other"** (static site).
5. Click **"Deploy"**.

---

## 🧪 Scientific Disclaimer & Safety Notice

This application is designed as an educational and computational platform. Model-derived outputs represent **predicted mathematical affinities** ($pIC_{50}$). They do not constitute clinical efficacy, biochemical proof of action, or safety in human beings. Every computationally prioritized candidate strictly requires rigorous wet-lab biochemical assays (e.g. TR-FRET kinase assays, cellular $IC_{50}$ assays in NSCLC cell lines) before advancement.
