# 4. Learning path: how we work through this together

Each session has a goal, something you read or do, and a concrete output
committed to the repo, so the understanding accumulates in Git and not just
in chat. Sessions build on each other, but 1–3 can be done before any
refactoring.

| # | Session | You read / do | We produce |
|---|---------|---------------|------------|
| 0 | **Environment** | Install R ≥ 4.4 and Rtools/Xcode. Then `install.packages(c("renv","mrgsolve","nlmixr2","tidyverse","survival","targets","tarchetypes","here","yaml","broom","marginaleffects","WeightIt","lmtp","gfoRmula"))` | `renv.lock`, `DESCRIPTION`. A smoke-test script that compiles each `.mod` |
| 1 | **The data-generating process** | The four `.mod` files and note 1 §1.3. Simulate *one* patient at 60 mg for each PK scenario and plot CP and tumour size | `analysis/00-model-walkthrough.qmd`: plots of PK1/PK2/PK3 accumulation, one tumour trajectory with and without drug, one DH2 dose history |
| 2 | **One script, end to end** | `ER1_DH1_1Dose.Rmd` with paths fixed. Run it at small N (200 IDs, 20 replicates) | Golden reference numbers (OR and p for Cavg1C and CavgTE). Verification of issues 1–2 in note 2 §2.3 (grid sensitivity of `INH` and of the AE draw) |
| 3 | **Why CavgTE lies** | Paper Fig. 2, plus the DAG in note 3 §3.1. Hand-compute CavgTE for an early versus a late event patient | A short explainer figure. You should be able to say *why* the slope is inverse under accumulation and positive under dose reduction |
| 4 | **Refactor I: functions** | Note 2 §2.4–2.5 | `R/exposure.R`, `R/endpoints.R`, `R/analyse_conventional.R` with unit tests. Golden numbers unchanged |
| 5 | **Refactor II: scenarios and pipeline** | `{targets}` walkthrough | `config/scenarios.yaml`, `_targets.R`. `tar_make()` reproduces the paper's Fig. 3 pattern at reduced replicates |
| 6 | **Estimands** | ICH E9(R1), plus note 3 §3.2 | Signed-off estimand table for the use case (we argue about intercurrent events here) |
| 7 | **Truth by counterfactual simulation** | Note 3 §3.3 | `R/truth.R`: true values of estimands A, B and C |
| 8 | **Causal estimators** | *What If* ch. 12–13 (IPW and standardisation), then ch. 19–21; the `lmtp` vignette | `R/analyse_causal.R`: g-computation, IPW, LMTP |
| 9 | **Model-based g-formula** | nlmixr2 fit of the joint model to one simulated trial | nlmixr2 → mrgsolve counterfactual pipeline, checked against the truth |
| 10 | **Simulation study and write-up** | ADEMP paper | `analysis/02-causal-use-case.qmd` with bias, coverage and power tables, and a recommendation for the dose question |

## Suggested first step

Sessions 0–2. Getting `ER1_DH1_1Dose.Rmd` to run unchanged except for
paths gives us the baseline everything else is checked against. It also
tests the two suspected model issues before we refactor around them.

## Conventions for our notes

* Explanations live in `docs/` (Markdown, rendered by GitHub).
* Runnable explanations live in `analysis/*.qmd` (Quarto), so the text and
  the numbers cannot drift apart.
* Every scenario is identified by its `id` from `config/scenarios.yaml`.
  When a note says "er2_dh2_500_200", that is reproducible with
  `tar_make(names = contains("er2_dh2_500_200"))`.
