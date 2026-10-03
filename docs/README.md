# Study notes: understanding, refactoring and extending this repository

These notes go with the code in `PSP-2025-0153/`. That code supports
Yin et al. (2025), *Evaluation and Mitigation of Time-Dependent Confounding
Effects in Conventional Exposure-Response Analyses for Oncology Drugs*,
CPT:PSP 14:2118–2127.

Read them in this order:

| # | File | What it gives you |
|---|------|-------------------|
| 1 | [`01-codebase-tour.md`](01-codebase-tour.md) | How the current code is built: the simulation pipeline, the models, the scripts and the order they depend on each other |
| 2 | [`02-issues-and-refactor-plan.md`](02-issues-and-refactor-plan.md) | Why it won't run out of the box, the technical problems we found, and a target layout that is configurable and reproducible |
| 3 | [`03-causal-roadmap.md`](03-causal-roadmap.md) | How to restate the question with estimands and causal inference, and a worked simulation use case to build in R (mrgsolve + nlmixr2 + tidyverse) |
| 4 | [`04-learning-path.md`](04-learning-path.md) | A session-by-session plan for working through all of this together |
| 5 | [`05-environment-setup.md`](05-environment-setup.md) | Where to run R (laptop, Codespaces, Hetzner, Claude cloud) and how to install everything |

The notes are opinionated on purpose. Where they disagree with the paper or
the code, they say why, and you are free to disagree back.
