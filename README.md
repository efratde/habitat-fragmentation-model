# Dispersal-Kernel Evolution under Habitat Fragmentation — an Agent-Based Model (NetLogo)

An individual-based (agent-based) simulation of how plant **seed-dispersal kernels** evolve when the landscape is **fragmented** into habitable and uninhabitable patches. Dispersal strategy is genetically encoded and free to change across generations, so the model lets selection reshape the dispersal kernel in response to the spatial structure of the habitat. Built in **NetLogo**, with parameter sweeps via **BehaviorSpace** and output analysis in **R**, documented following the **ODD protocol**.

> **Note on provenance — please read.** This is **my own PhD agent-based model** of dispersal-kernel evolution **under habitat fragmentation** — one of a family of dispersal-evolution models I designed and built during my doctorate. It is an exploratory research model of mine and was **not** the basis of any specific publication.
>
> **This is NOT the model used in Greenbaum, Dener & Giladi 2022** (*J. R. Soc. Interface* — "Limits to the evolution of dispersal kernels under rapid fragmentation"). That paper used a **different model, developed by G. Greenbaum**, and is listed on my CV under publications. The two models share the *word* "fragmentation" and a broad research theme, but they are **separate pieces of software**: the code in this repository is mine and is unrelated to the code behind that paper. Please do not let the shared word "fragmentation" lead you to conflate them.

---

## What the model does

Plants disperse seeds according to a *dispersal kernel* — the probability distribution of dispersal distances. This model simulates a population of plants whose dispersal strategy is **genetically encoded** (heritable alleles), living on a spatial landscape that can be **fragmented** into habitable and uninhabitable patches. Seeds that land on uninhabitable patches fail to establish, so the geometry of fragmentation imposes selection on the dispersal kernel. The simulation runs dispersal, pollination, establishment, competition, and reproduction across many generations, tracking how the dispersal kernel — and the population's **spatial genetic structure** and **inbreeding load** — responds to the fragmented habitat.

A `fragmentation` chooser controls the landscape regime:

- **`no`** — a continuous, fully habitable landscape (no fragmentation);
- **`black uninhabitable`** — dark patches are non-habitat, breaking the landscape into fragments;
- **`white uninhabitable`** — light patches are non-habitat, the complementary fragmentation pattern.

Sweeping parameters (fragmentation regime, landscape autocorrelation, dispersal distance, gamete number, neighbourhood size, genetic architecture, mutation rate) maps how the evolving dispersal strategy and the population's genetic structure co-respond to habitat fragmentation.

## My role

- **Designed and implemented the model** in NetLogo — population dynamics, genetic encoding of dispersal strategies, the fragmented-landscape generator, and the dispersal / pollination / establishment / competition / reproduction rules.
- **Designed and ran the simulation experiments** (BehaviorSpace parameter sweeps over fragmentation regimes and dispersal parameters).
- **Analyzed the output in R** and produced summary figures.
- Documented the model following the **ODD protocol** for reproducibility.

## Repository structure

```
habitat-fragmentation-model/
├── model/        NetLogo model (.nlogo) — the full desktop simulation, with the
│                 BehaviorSpace parameter-sweep experiments embedded inside it;
│                 also a Lite *-web.nlogo* variant that runs in NetLogo Web (browser)
└── .gitignore    excludes simulation output, raw data, and bulky files
```

## Data policy

This repository contains **code, not data**. The simulation **generates** its own data — there is nothing private here. The model reads **no external files**: the only file it touches is the run output that **BehaviorSpace writes itself** (one CSV per run). That raw output (CSV/XLSX/RData) is regenerable and is excluded via `.gitignore` to keep the repo lean.

## Running the model

This is a **NetLogo 6.x desktop** model. It relies on the `rnd`, `stats`, and `profiler` extensions, which **NetLogo Web does not support** — so the full model runs in the desktop application, **not** in the browser. A reduced **Lite** variant, `model/habitat-fragmentation-dispersal-model-web.nlogo`, strips those extensions (the one weighted-random draw is reimplemented in pure NetLogo, and the `stats`/`profiler`-driven analysis widgets are removed) so the core dispersal-evolution dynamics can run directly in **NetLogo Web** — see *Try it in the browser* below.

1. Open `model/habitat-fragmentation-dispersal-model.nlogo` in **[NetLogo 6.x](https://ccl.northwestern.edu/netlogo/)**, with the `rnd`, `stats`, and `profiler` extensions available.
2. Choose a landscape regime with the **fragmentation** chooser (`no`, `black uninhabitable`, or `white uninhabitable`).
3. Click **setup**, then **go** to run the simulation interactively.
4. For experiments: run the embedded BehaviorSpace experiment (`Tools → BehaviorSpace`), then open the output in R/RStudio to reproduce the figures.

### Try it in the browser (Lite)

No install needed: open the Lite variant in **[NetLogo Web](https://www.netlogoweb.org/launch)** (choose *Upload a Model* → `model/habitat-fragmentation-dispersal-model-web.nlogo`), then click **setup** and **go**. The Lite variant keeps the dispersal / pollination / establishment / reproduction dynamics intact but drops the `stats`/`profiler`-driven analysis widgets and the BehaviorSpace sweeps, which live only in the full desktop model above.

---

*Part of the research-software portfolio of Dr. Efrat Dener — plant ecologist & quantitative/computational researcher.*
