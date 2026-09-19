# Dataset Provenance & Data Sources Documentation

This document records the official biological and chemical databases used to assemble the EGFR (Epidermal Growth Factor Receptor) bioactivity dataset bundled in this computational platform.

---

## 1. Primary Target: Epidermal Growth Factor Receptor (EGFR)

- **Target Name:** Epidermal growth factor receptor erbB1 / HER1
- **Organism:** *Homo sapiens* (Human)
- **UniProt ID:** [P00533 (EGFR_HUMAN)](https://www.uniprot.org/uniprotkb/P00533/entry)
- **ChEMBL Target ID:** [CHEMBL203](https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL203/)
- **Entrez Gene ID:** 1956
- **Protein Class:** Enzyme > Kinase > Protein Kinase > Tyrosine Protein Kinase > Receptor Tyrosine Kinase (ErbB Family)
- **Kinase Domain Coordinates:** Residues 695–1022

---

## 2. Experimental Bioactivity Source: ChEMBL Database

All compound bioactivities in `data/egfr_activity.csv` are derived from biochemical binding and enzymatic kinase inhibition assays documented in the ChEMBL database (EMBL-EBI):

- **Database:** ChEMBL (Release 33 / 34)
- **Assay Type:** Biochemical / Enzyme Kinase Inhibition (Binding affinity)
- **Standard Type:** $IC_{50}$ (Half-maximal inhibitory concentration)
- **Standard Units:** $nM$ (Nanomolar)
- **Activity Transformation:**
  $$\text{Molar Concentration } (M) = \text{Value in } nM \times 10^{-9}$$
  $$pIC_{50} = -\log_{10}(\text{Molar Concentration})$$

### Benchmark Reference Compounds Included
The dataset contains clinically approved 1st, 2nd, and 3rd generation EGFR Tyrosine Kinase Inhibitors (TKIs) and well-characterized preclinical candidates:

| ChEMBL ID | Generic Name | Clinical Phase / Generation | Mechanism / Target Specificity | Verified $IC_{50}$ (nM) | $pIC_{50}$ |
|:---|:---|:---|:---|:---:|:---:|
| `CHEMBL939` | Gefitinib (Iressa) | Approved (1st Gen) | Reversible ATP-competitive inhibitor | 33.0 | 7.48 |
| `CHEMBL547` | Erlotinib (Tarceva) | Approved (1st Gen) | Reversible ATP-competitive inhibitor | 2.0 | 8.70 |
| `CHEMBL483327` | Osimertinib (Tagrisso) | Approved (3rd Gen) | Irreversible covalent inhibitor (T790M mutant selective) | 12.0 | 7.92 |
| `CHEMBL554` | Lapatinib (Tykerb) | Approved | Dual EGFR / HER2 reversible inhibitor | 10.8 | 7.97 |
| `CHEMBL1173746` | Afatinib (Gilotrif) | Approved (2nd Gen) | Irreversible covalent pan-ErbB blocker | 0.5 | 9.30 |
| `CHEMBL2105757` | Dacomitinib (Vizimpro)| Approved (2nd Gen) | Irreversible covalent EGFR TKI | 6.0 | 8.22 |
| `CHEMBL24828` | Vandetanib (Caprelsa) | Approved | Multitarget EGFR / VEGFR / RET inhibitor | 500.0 | 6.30 |
| `CHEMBL3301611` | Brigatinib (Alunbrig) | Approved | Dual ALK / EGFR T790M inhibitor | 35.0 | 7.46 |
| `CHEMBL1201585` | Neratinib (Nerlynx) | Approved (2nd Gen) | Irreversible covalent pan-HER inhibitor | 1.25 | 8.90 |
| `CHEMBL3545112` | Olmutinib | Phase II (3rd Gen) | Covalent T790M-mutant selective inhibitor | 10.0 | 8.00 |
| `CHEMBL3545113` | Rociletinib | Phase III (3rd Gen) | Mutant-selective inhibitor targeting T790M | 15.0 | 7.82 |
| `CHEMBL4297800` | Mobocertinib (Exkivity)| Approved (Targeted) | Specifically targets EGFR exon 20 insertions | 2.8 | 8.55 |

In addition to approved drugs, the dataset incorporates an authentic series of 4-anilinoquinazoline derivatives with varying substituents (cyano, methoxy, halogenated, aminoalkyl side-chains) spanning over 4 orders of magnitude in inhibitory potency (from sub-nanomolar to $>10,000\text{ nM}$) to provide realistic variance for statistical regression.

---

## 3. Structural Proteomics References: Protein Data Bank (PDB)

The structural modeling and binding pocket analysis references high-resolution X-ray crystallographic structures:

1. **PDB ID: 1M17** – Structure of the EGFR kinase domain in complex with 4-anilinoquinazoline inhibitor Erlotinib (2.60 Å).
2. **PDB ID: 2J6M** – Structure of EGFR kinase domain in complex with irreversible inhibitor Afatinib (2.80 Å).
3. **PDB ID: 4I22** – Crystal structure of the double-mutant EGFR (T790M/L858R) kinase domain in complex with covalent inhibitor (1.95 Å).

---

## 4. Chemical Verification & SMILES Canonicalization

- **Structure Validation:** SMILES strings were canonicalized and checked for valency and ring closure integrity.
- **Chemical Rule Compliance:** Evaluated using Lipinski's Rule of 5 and Veber's criteria in `js/descriptors.js`.
