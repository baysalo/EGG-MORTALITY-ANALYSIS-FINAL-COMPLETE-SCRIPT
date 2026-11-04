# Egg Mortality Statistical Analysis Pipeline (17-Group Experimental Design)

An R-based analytical workflow for evaluating egg mortality across 16 treatment groups and 1 untreated control group. 

This pipeline processes grouped-binomial response data using **bias-reduced binomial logistic regression**, performs an **omnibus replicate-level permutation test**, conducts **exact pairwise permutation post-hoc tests with Holm multiplicity adjustment**, and generates a **compact letter display (CLD)** alongside automated validation checks.

---

## 📌 Experimental Design

- **Treatment Groups:** 17 groups total (16 treatments + 1 untreated control designated as `Treatment 17`).
- **Replication:** 8 independent Petri-dish replicates per treatment group ($N = 136$ total experimental units).
- **Sample Size:** 10 eggs per Petri dish ($N = 1,360$ total eggs).
- **Primary Endpoint:** Mortality percentage (`%mortality`), converted internally to discrete `Dead` and `Alive` counts.

---

## 🔬 Statistical Methodology

1. **Data Pre-processing & Validation:**
   - Enforces structural data checks (17 groups, 8 replicates/group, 10 eggs/replicate).
   - Converts mortality percentage into exact discrete integer counts of `Dead` and `Alive` eggs.
   - Validates mathematical identity ($Dead + Alive = 10$) and checks for duplicate identifiers.

2. **Bias-Reduced Binomial Logistic Regression:**
   - Fits a generalized linear model (GLM) with a logit link function using **Firth’s mean bias-reduction** method via `brglm2::brglmFit` (`AS_mean` type).
   - Prevents parameter inflation caused by complete or quasi-complete separation (e.g., 0% or 100% mortality).

3. **Omnibus Replicate-Level Permutation Test:**
   - Evaluates overall treatment effect while preserving whole Petri-dish units (136 independent observations).
   - Computes an empirical $p$-value over $B = 100,000$ Monte Carlo permutations using a binomial Likelihood Ratio (LR) test statistic.

4. **Exact Pairwise Post-Hoc Permutation Tests:**
   - Evaluates all $\binom{17}{2} = 136$ pairwise group comparisons.
   - For each pair, evaluates all exact permutations ($\binom{16}{8} = 12,870$ allocations per pair) across $1,750,320$ total exact iterations.
   - Adjusts raw exact $p$-values using the **Holm-Bonferroni** step-down method ($\alpha = 0.05$).

5. **Compact Letter Display & Verification:**
   - Grouping via `multcompView::multcompLetters` based on Holm-adjusted $p$-values.
   - Includes an automated consistency check verifying that groups sharing a letter do not significantly differ, and non-overlapping groups differ significantly.

---

## 🛠️ Requirements & Installation

### R Dependencies

The script automatically detects missing CRAN packages and installs them upon execution. The required libraries are:

```R
install.packages(c("brglm2", "readxl", "dplyr", "multcompView"), dependencies = TRUE)
```

---

## 📂 Data Input Format

The script prompts the user to select an Excel spreadsheet interactively (`file.choose()`). 

The Excel file **must** contain a sheet named **`yumurta`** with the following column headers:

| Column Header | Type | Description |
| :--- | :--- | :--- |
| `Treat` | Numeric / Character | Group identifier (1 to 17; `17` = Untreated Control) |
| `Replication` | Integer | Replicate identifier (1 to 8 per group) |
| `%mortality` | Numeric | Percent egg mortality (e.g., `0`, `10`, `100`) |

---

## 🚀 Execution Guide

1. Clone or download this repository:
   ```bash
   git clone https://github.com/your-username/egg-mortality-analysis.git
   cd egg-mortality-analysis
   ```

2. Open R or RStudio and source the script:
   ```R
   source("egg_mortality_analysis.R")
   ```

3. When prompted by the file selection dialog, select your `.xlsx` input file.

---

## 📊 Output Files Generated

Upon completion, the script automatically exports five formatted CSV files into the working directory:

| Filename | Contents |
| :--- | :--- |
| `Egg_17groups_descriptive_results.csv` | Group-level summary (Dead, Alive, Total, Observed %, Predicted %). |
| `Egg_17groups_omnibus_permutation.csv` | LR statistic, Monte Carlo SE, and permutation $p$-value. |
| `Egg_17groups_exact_pairwise_Holm.csv` | Full pairwise table (136 pairs) with unadjusted and Holm-adjusted $p$-values. |
| `Egg_17groups_significant_pairs.csv` | Filtered subset of statistically significant pairwise comparisons ($\alpha < 0.05$). |
| `Egg_17groups_Holm_letters.csv` | Final compact letter display table for publication/reporting. |
| `Egg_17groups_letter_check.csv` | Audit trail for letter assignment validity. |

---

## 📜 Session Info & Reproducibility

- **Random Seed:** Set to `20261029` for exact reproducibility of the omnibus permutation test.
- **Monte Carlo Iterations:** 100,000 permutations.
- Detailed R session information is printed to the console at script completion.

  To cite: Mehmet Kececi and Omur Baysal, 2026 
