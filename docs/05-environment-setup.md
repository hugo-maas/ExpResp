# 5. Setting up the R environment

## 5.1 Where should the work happen?

Code and compute are separate questions. **GitHub is where the code lives**:
the history, the branches and our notes. It is not where R runs, unless you
use GitHub Codespaces. Four places can run R for this project, and you will
probably use more than one:

| Where | Good for | Cost | Effort | Verdict |
|---|---|---|---|---|
| **Your own computer** (Positron or RStudio) | Learning sessions 0–8: interactive plots, stepping through code, small N | free | low–medium | **Start here** |
| **GitHub Codespaces** (browser VS Code in a container) | Working from any machine with identical setup; quick experiments | free tier currently about 120 core-hours/month on personal accounts, then pay per hour | low once a dev container exists | Nice to have |
| **Hetzner Cloud VM** | Replicate sweeps (sessions 9–10): hundreds of replicates × scenarios on 16–48 cores; also running the *original* hourly-grid scripts, which need a lot of RAM | billed hourly, roughly tens of € per month for a 16-vCPU machine if left running; check current prices | medium | **Add when simulations get big** |
| **The Claude Code cloud environment** (this session) | Letting Claude run, test and debug the R code itself rather than only reading it | included | one-time settings change | **Enable now** (§5.5) |

My recommendation:

1. Install R locally now. You need to see plots and poke at objects to
   build understanding, and a laptop is the best place for that.
2. At the same time, enable R in this Claude environment (§5.5), so I can
   verify the code changes I make before pushing.
3. Rent a Hetzner machine only when the replicate studies start, and delete
   it afterwards. The `{targets}` pipeline caches results in `_targets/`,
   so runs can be resumed and their results copied home.
4. Make all of these identical by committing `renv.lock`, the exact package
   versions, to the repo. A dev container (`.devcontainer/`) is optional
   and becomes useful once you use Codespaces or Hetzner.

### How much compute is needed?

The original scripts simulate an **hourly grid × 600 days × 1000
subjects**. That is about 14.4 M rows and about 25 columns, so roughly
**3 GB of RAM per simulation** before any copies, and they repeat it for
each replicate. 16 GB of RAM is the practical minimum for running them
unchanged, and 32 GB is comfortable. After the refactor (daily grid,
note 2 §2.3 point 6), a 1000-subject trial should take seconds and fit
easily on a laptop. Many replicates are then a **CPU** problem, which is
where a many-core VM pays off.

## 5.2 Local installation, step by step

### 1. R itself

Use the current R release (4.5.x or newer; the paper used 4.4.0) from
<https://cloud.r-project.org>.

### 2. A C/C++ toolchain

This is required: **mrgsolve compiles every model to C++**, and nlmixr2
and rxode2 compile too.

| OS | What to install | Check |
|---|---|---|
| Windows | **Rtools** matching your R version (e.g. Rtools45 for R 4.5) from CRAN. Do *not* hard-code the PATH as the original scripts do; the installer handles it | `pkgbuild::has_build_tools(debug = TRUE)` |
| macOS | `xcode-select --install` in Terminal, plus the gfortran build listed on CRAN's macOS tools page | same |
| Linux (Ubuntu) | `sudo apt install build-essential gfortran libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev` | same |

### 3. An IDE

* **Positron** (Posit's newer IDE, VS Code based, good Git and Quarto
  support) or **RStudio Desktop**. Either works.
* Install **Quarto** (<https://quarto.org>) for the `analysis/*.qmd`
  documents.

### 4. Git and the repository

```bash
git clone https://github.com/hugo-maas/ExpResp.git
cd ExpResp
git switch claude/exposure-response-refactor-i63khh
```

Authenticate once with `gh auth login` (GitHub CLI) or an SSH key. Then open
the folder as a project in Positron or RStudio, so the working directory is
always the repo root. This is the habit that replaces `setwd()`.

### 5. Packages, pinned with renv

In R, from the repo root:

```r
install.packages("renv")
renv::init(bare = TRUE)   # project-private library, records versions

# Fast binary installs on Linux / Windows / macOS:
options(repos = c(CRAN = "https://p3m.dev/cran/latest"))

renv::install(c(
  # core
  "tidyverse", "here", "yaml", "glue",
  # simulation & estimation
  "mrgsolve", "nlmixr2", "rxode2",
  # conventional E–R as in the paper
  "survival", "survminer", "broom", "zoo", "patchwork", "openxlsx",
  # causal inference
  "marginaleffects", "WeightIt", "lmtp", "gfoRmula", "SuperLearner",
  # workflow
  "targets", "tarchetypes", "crew", "quarto", "testthat", "devtools"
))
renv::snapshot()          # writes renv.lock → commit it
```

* **`brms`** (used for the paper's Emax logistic fits) needs Stan through
  `cmdstanr` or `rstan`. Leave it out until we need it; it is the
  heaviest dependency.
* **nlmixr2** is large. On Linux, install binaries from p3m (as above) or
  r2u (§5.4). Compiling everything from source can take 30–60 minutes.

### 6. Smoke test

These lines confirm that the toolchain and the key packages work:

```r
library(mrgsolve)
mod <- mread("KP_2com", project = "PSP-2025-0153", file = "KP_2com.mod")
out <- mrgsim(mod, ev(amt = 60, ii = 24, addl = 13), end = 24 * 14, delta = 1)
plot(out, CP ~ time)            # should show accumulation toward steady state

library(nlmixr2)
one.cmt <- function() {
  ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
  model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
          linCmt() ~ add(add.sd) })
}
fit <- nlmixr2(one.cmt, theo_sd, est = "focei", control = list(print = 0))
fit$objf                        # a number = estimation works
```

If both run, session 0 is done.

## 5.3 GitHub Codespaces (optional)

When we add `.devcontainer/devcontainer.json` based on a **rocker** image
(e.g. `rocker/r-ver` with the toolchain, Quarto and the R extension), you
open the repo on GitHub and choose **Code → Codespaces → Create**. You get
VS Code in the browser with R ready. Choose a 4-core machine or larger.
Stop the codespace when you are done, because it bills by the hour above
the free quota. The same container definition also works locally with
Docker and on Hetzner, which is the main reason to create it.

## 5.4 Hetzner Cloud VM (for the big runs)

1. **Create the server.** In the Hetzner Cloud console:
   * Choose Ubuntu 24.04.
   * Choose a **dedicated-vCPU** type (CCX line). Simulations are
     CPU-bound, and shared vCPUs throttle. 16 vCPU and 64 GB is a good
     size.
   * Add your **SSH key**. Do not use password login.
   * Attach a **firewall** allowing only port 22.
2. **Install R from binaries** with r2u, which turns CRAN into Ubuntu
   packages so installs take minutes:
   ```bash
   # as root on the server
   curl -fsSL https://raw.githubusercontent.com/eddelbuettel/r2u/master/inst/scripts/add_cranapt_noble.sh | bash
   apt install -y r-cran-tidyverse r-cran-mrgsolve r-cran-nlmixr2 r-cran-targets \
                  r-cran-crew r-cran-lmtp r-cran-weightit r-cran-survival git
   ```
   Alternatively, run `renv::restore()` with the p3m repo, so that the
   versions match your laptop exactly.
3. **Connect without exposing anything.**
   * Positron or VS Code **Remote-SSH** to the server, or
   * RStudio Server through an SSH tunnel:
     `ssh -L 8787:localhost:8787 you@server`, then open
     <http://localhost:8787>. Never open port 8787 publicly.
4. **Run long jobs so they survive disconnects.** Use `tmux`, then
   `Rscript -e 'targets::tar_make()'`. `{crew}` uses all cores.
5. **Stop paying.** Hetzner bills a server while it exists, **even when
   powered off**. When done, copy `_targets/` and the outputs back
   (`rsync`), take a snapshot if you want to resume later, and **delete**
   the server.

Data: the repo contains only simulated data and published PK parameter
estimates, so an EU cloud VM raises no special data-protection concerns. If
you ever bring in real trial data, that changes.

## 5.5 Enabling R in this Claude cloud environment

Right now I can read your code but not run it. This environment's network
policy blocks the R package servers. I tested `cloud.r-project.org`,
`p3m.dev`, `packagemanager.posit.co` and `r2u.stat.illinois.edu`, and all
were denied. Ubuntu's own R is reachable but old (4.3.3) and comes without
CRAN packages.

To fix it, open the **cloud environment menu in the session title bar →
Edit**:

1. Under **Network access**, choose *Custom*, keep the default package
   managers, and add these allowed domains:
   * `cloud.r-project.org`
   * `r2u.stat.illinois.edu`
   * `p3m.dev`
   * `keyserver.ubuntu.com`
   * `raw.githubusercontent.com`

   Docs: <https://code.claude.com/docs/en/cloud-environments#network-access>
2. Under **Setup script**, paste:
   ```bash
   #!/bin/bash
   set -euo pipefail
   # R + CRAN packages as Ubuntu binaries (r2u)
   curl -fsSL https://raw.githubusercontent.com/eddelbuettel/r2u/master/inst/scripts/add_cranapt_noble.sh | bash
   apt-get install -y --no-install-recommends \
     r-cran-tidyverse r-cran-mrgsolve r-cran-nlmixr2 r-cran-rxode2 \
     r-cran-survival r-cran-survminer r-cran-broom r-cran-zoo r-cran-patchwork \
     r-cran-openxlsx r-cran-here r-cran-yaml r-cran-targets r-cran-tarchetypes \
     r-cran-crew r-cran-testthat r-cran-devtools r-cran-marginaleffects \
     r-cran-weightit r-cran-lmtp r-cran-gformula r-cran-superlearner
   ```
   This script has **not been tested yet** because the domains are
   blocked. New sessions run it; if it fails, the session shows which
   command failed and I'll adjust it.

These cloud sessions have about 4 cores and 15 GB of RAM: enough for
development and small-N checks, not for the full replicate sweeps.
