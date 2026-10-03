# 2. Why it doesn't run out of the box, and how to fix it

## 2.1 My overall opinion

The *science* is sound and useful. It makes a clear, simulation-backed
argument that metrics like CavgTE are contaminated by the outcome process,
so they are contaminated by design.

The *code* is a researcher's lab notebook, not a tool. It records what one
person did on one machine, interactively, while changing parameters by hand
between runs. That is normal for a paper's supplement, but it means:

* **Nothing runs as-is.** Every script `setwd()`s into
  `C:/Users/Xuefen.Yin/OneDrive - FDA/...`, and several read intermediate
  files that only exist there.
* **The "parameters" are edits.** Scenarios are defined by changing numbers
  in the middle of chunks and renaming output folders by hand. The repo
  therefore cannot tell you which settings produced which figure.
* **It is about 20,000 lines with roughly 90 % duplication.** The same
  helper functions and pipeline appear in 19 files with drifting
  differences.

The good news is that the underlying pipeline is small. Once deduplicated
it is perhaps **5 or 6 functions and 3 model files**, plus a scenario table.
That structure is what makes "different parameter settings out of the box"
easy.

## 2.2 Blocking problems (must fix to run at all)

| # | Problem | Where | Fix |
|---|---------|-------|-----|
| B1 | 324 `setwd()` calls and 762 lines with absolute `C:/Users/...` paths (`setwd`, `readRDS`, `ggsave`) | every Rmd | `here::here()` and one `output_dir` derived from the scenario ID |
| B2 | `dosing record for study 309.csv` (real clinical dose records) not provided | all 4 `ER1_DH3_*` | Replace DH3 with a **synthetic empirical-dosing generator** (for example a discrete-time Markov chain over dose levels and interruption, calibrated to the pattern in paper Fig. S2) |
| B3 | `PKPD_tumor_DH2_AE.mod` referenced but absent | `ER1_DH2_2Dose_ORR_TTE.Rmd` | Probably the same as `…_LDR.mod`; confirm against the paper and consolidate |
| B4 | Four model files and the "true ER" CSVs exist only inside `True ER relationship.zip` | true-ER and ER2 PFS scripts | Unzip into `models/` and `reference/`, or better, regenerate the truth in code (§2.5) |
| B5 | Cross-script dependencies on intermediate files (`df_exp.rds`, `ORR_df_ER_IDQ_*.csv`) with no documented run order | see the graph in note 1 §1.5 | A pipeline tool (`{targets}`) that encodes the dependency graph |
| B6 | Windows-specific Rtools PATH manipulation | setup chunks | Remove it. Toolchains are set up per OS, not in analysis code |
| B7 | `rm(list = ls())` at the top of notebooks | setup chunks | Remove it. Knit in a clean session instead |

## 2.3 Correctness and robustness issues

These were found by reading, not by running (R is not yet installed in
this environment), so treat them as **hypotheses to verify in session 2**.

1. **The drug effect in the ODE is evaluated in `$ERROR`.**
   `INH = Emax·CPⁿ/(CPⁿ+EC50ⁿ)` is computed in `$ERROR` but used in
   `$ODE`. mrgsolve promotes the declaration to a global, so the ODE sees
   the value from the *last output record*: a piecewise-constant, lagged
   drug effect whose accuracy depends on `tgrid`. With the hourly grid used
   here the error is probably small, but results change if you coarsen the
   grid. **Fix:** compute `INH` inside `$ODE` from `CENT`.
2. **The AE model draws random numbers on every output record.** The
   `R::rbinom()` call in `$ERROR` runs every hour, but only the draw at
   `TIME %% 24 == 0` is acted on. The random stream, and so the dose
   history, therefore depends on the output grid. Changing `tgrid` changes
   which patients get dose reductions, even with the same seed. **Fix:**
   draw the AE once per day at the dosing decision, or simulate the
   dose-modification process outside the ODE in R, day by day.
3. **Hand-set parameters disagree with the folder they are saved to.** For
   example `ER2_PFS_DH2_AE_LDR_1D_1000ID.Rmd` sets doses `1 / 0.6 / 0.4`
   but saves to `…_500&200&100_LDR/`. `ER2_ORR_DH2…` sets `500/200/100`
   labels and then reloads the `60&40&20` file. Some figures may be
   mislabelled, and we cannot know which from the code alone.
4. **ER1 uses a 10 mg dose while the paper text says 60 mg.** Under linear
   PK this only rescales the x-axis and does not change the conclusions,
   but it is confusing.
5. **Per-ID `for (i in 1:1000)` loops assume `ID == row index`.** For
   example `Analysis.PFS$PFS[i]` is paired with `df.exp$ID == i`, which
   breaks silently if IDs are reordered or subset. **Fix:** use a `left_join`
   and `filter(DAY <= event_day)`. This is also about 100× faster.
6. **Memory.** An hourly grid × 600 days × 1000 IDs ≈ 14.4 M rows × about
   25 columns, all held in RAM, and repeated per replicate. Exposure
   metrics need only daily AUC increments, plus dense sampling in cycle 1
   if you want Cmax. **Fix:** a daily `tgrid`, with mrgsolve computing
   cumulative AUC in a compartment and `Cmax`/`Ctrough` via captured
   variables or a cycle-1 dense grid.
7. **Endpoint definitions are simplified.** PFS is defined against baseline
   (T ≥ 1.2·T0) rather than the RECIST nadir (+20 % and +5 mm). Tumours
   are "assessed" daily rather than every 6–8 weeks, and ORR needs no
   confirmation. This is fine for the paper's purpose, but a realistic use
   case should add an **assessment schedule**: it creates interval
   censoring, which interacts with time-dependent exposure.
8. **No baseline confounding exists in the data-generating model.** PK
   parameters are independent of tumour parameters. That is deliberate, to
   isolate the time-dependent mechanism, but it means a *static* cycle-1
   exposure is automatically unconfounded here. In real data (CL
   correlated with tumour burden, albumin, cachexia, as with nivolumab or
   pembrolizumab) it is not. This is the main gap for the causal work
   (note 3).

## 2.4 The target layout

```
ExpResp/
├── DESCRIPTION              # makes it an R "project package": deps + devtools::load_all()
├── renv.lock                # exact package versions
├── _targets.R               # the pipeline (dependency graph, caching, parallel)
├── config/
│   └── scenarios.yaml       # ONE row/entry per scenario = every knob in one place
├── models/                  # .mod / .cpp for mrgsolve (one canonical PK-TGI-AE model)
├── data-raw/                # PK parameter CSV, synthetic DH3 generator inputs
├── R/
│   ├── simulate.R           # sim_population(), sim_trial()  (mrgsolve)
│   ├── dosing.R             # dosing_dh1(), dosing_dh2_rule(), dosing_dh3_markov()
│   ├── exposure.R           # derive_exposure(sim, metrics = c("Cavg1C","CavgSS","CavgTE",...))
│   ├── endpoints.R          # derive_orr(), derive_pfs(), assessment schedule, truncation
│   ├── er1_outcome.R        # Weibull TTE + censoring, independent of drug
│   ├── analyse_conventional.R  # logistic, Emax, KM-by-quantile, Cox -> tidy tibble
│   ├── analyse_causal.R     # g-computation, IPW/MSM, LMTP (note 3)
│   ├── truth.R              # compute the true estimand by counterfactual simulation
│   └── plots.R              # one set of plotting functions
├── analysis/                # Quarto documents that *read* targets results and explain them
│   ├── 01-reproduce-paper.qmd
│   └── 02-causal-use-case.qmd
├── tests/testthat/          # unit tests + "golden" checks against paper figures
└── PSP-2025-0153/           # the original code, left untouched for reference
```

### What a scenario looks like

```yaml
# config/scenarios.yaml
defaults:
  n_per_arm: 500
  days: 600
  n_rep: 200
  seed: 20871
  pk_scenario: PK2          # PK1 | PK2 | PK3 -> KA/CL/V multipliers live in one table
  dosing: DH1               # DH1 | DH2 | DH3
  interval_h: 24
  er: ER2                   # ER1 (Weibull, null) | ER2 (tumour model)
  endpoints: [ORR, PFS]
  exposure_metrics: [Cavg1C, CavgSS, CavgTE, mCavgTE]

scenarios:
  - id: er1_dh1_pk3q2w_eo1
    er: ER1
    pk_scenario: PK3
    interval_h: 336
    arms: [{start_dose: 60}]
    weibull: {shape: 2, scale: 110}
    dropout: 0.25

  - id: er2_dh2_500_200
    dosing: DH2
    arms: [{start_dose: 500, reductions: [200, 100]},
           {start_dose: 200, reductions: [100, 60]}]
```

Running the code is then:

```r
# once
renv::restore()

# everything (cached; only changed scenarios rerun)
targets::tar_make()

# or one scenario interactively, with overrides
devtools::load_all()
sc  <- read_scenario("er2_dh2_500_200", n_rep = 20, n_per_arm = 200)
sim <- sim_trial(sc)
exp <- derive_exposure(sim, sc$exposure_metrics)
res <- analyse_conventional(exp, endpoints = "ORR")
```

`{tarchetypes}::tar_map()` turns each YAML entry into its own branch of the
pipeline, and `{crew}` runs replicates in parallel. Output folders are
named from `scenario$id` automatically, so the label-and-content mismatch
described in §2.3 point 3 cannot happen again.

## 2.5 Refactoring strategy: preserve behaviour first, then improve

1. **Freeze the reference.** Keep `PSP-2025-0153/` untouched.
2. **Get one script running end to end with only path fixes.** Use
   `ER1_DH1_1Dose.Rmd`, since it has no hidden inputs. Record its key
   outputs (Type I error and OR for CavgTE/Cavg1C) as a *golden file*.
3. **Extract functions one stage at a time.** Order: exposure, endpoints,
   analysis, simulation. After each step, check that the golden numbers
   are unchanged with the same seeds. Only *after* that, apply the
   behaviour-changing fixes (INH in `$ODE`, daily AE draw, daily grid), and
   document how much the numbers move.
4. **Generalise.** Merge the four model files into one parameterised model
   with switches (`AE_ON`, `PD_ON`, `EMAX`), and unify DH1/DH2/DH3 as
   dosing generators.
5. **Acceptance tests.** Reproduce the *direction* of the paper's key
   findings (Fig. 2a/b: inverse slope under accumulation; Fig. 2c/d:
   positive slope under DH3; Fig. 3: static metrics keep Type I error ≤
   10 %). Use these as `testthat` checks on small replicate counts.

## 2.6 Tooling notes for R

* **mrgsolve** stays the simulation engine. Use `mread()` from `models/`
  with `project = here("models")`, plus `obsonly()`, `outvars` and a daily
  `tgrid`.
* **nlmixr2** (the current name of nlmixr) is for the *estimation* side:
  fitting the PK/tumour model to a simulated trial as an analyst would
  (note 3). Models are written once in mrgsolve syntax and once as an
  nlmixr2/rxode2 function. A unit test should simulate both with the same
  parameters and check that the profiles agree.
* **tidyverse** throughout, with the base-R `for` loops and `merge` calls
  replaced.
* **survival**, plus **broom** for tidy model output, instead of the
  hand-written `univ_con` / `univ_cat` / `univ_cat_T` extractors.
* **brms** for the Emax logistic model is slow inside replicate loops.
  Consider `nls`/`gnm`-style maximum likelihood (or `glm` with a fixed
  Hill term) for replicates, and keep brms for showcase fits.
