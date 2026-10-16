# MMDP

**Integrating complementary datasets to construct microbiome–metabolome–disease pathways**

`MMDP` implements the method of Fang, Ye and Wang (Bioinformatics, accepted). It studies how gut microbes,
metabolites and a disease outcome are linked when the data are split across two studies:

- a **target** dataset: a large clinical cohort with microbiome profiles and the outcome, but **no metabolomics**;
- an **external** dataset: a smaller study with paired microbiome and metabolite profiles, but without the outcome
  (or with too few subjects to study it).

The package predicts each metabolite in the target cohort from the microbiome and tests its association with the
outcome while adjusting for microbial confounders. Predicting metabolites from the microbiome across datasets
has worked well in earlier studies (Mallick et al. 2019). The package also sorts microbes into functional
groups. To avoid
results that depend on an arbitrary choice of reference taxon in the additive log-ratio (ALR) transform, every
microbe is used as the reference in turn and the results are combined.

---

## Model

For subject $i$ in the target study, let $\mathbf{w}_i$ be the abundances of $p$ taxa, $m_i$ the unmeasured
metabolite and $y_i$ the outcome. The structural equation model is

$$
m_i = (\log \mathbf{w}_i)^\top \boldsymbol\gamma^0 + \delta_i, \qquad \mathbf{1}_p^\top\boldsymbol\gamma^0 = 0,
$$

$$
y_i = m_i\ \theta^* + (\log \mathbf{w}_i)^\top \boldsymbol\beta^0 + \varepsilon_i, \qquad \mathbf{1}_p^\top\boldsymbol\beta^0 = 0,
$$

where $\theta^*$ is the direct effect of the metabolite on the outcome. The errors $\delta_i$ and $\varepsilon_i$
may be correlated, which allows for unmeasured confounding such as diet. The zero-sum constraints follow from the
log-contrast model for compositional data (Aitchison 1982).

The coefficient pair $(\gamma^0_l, \beta^0_l)$ puts each microbe $l$ in one of four functional groups:

| Group | $\gamma^0_l$ | $\beta^0_l$ | Interpretation | `stability_selection()` output |
|---|---|---|---|---|
| $\mathcal G_1$ | $\neq 0$ | $\neq 0$ | associated with both the metabolite and the outcome | `$both` |
| $\mathcal G_2$ | $\neq 0$ | $= 0$ | affects the outcome only **through** the metabolite | `$stage1_only` |
| $\mathcal G_3$ | $= 0$ | $\neq 0$ | affects the outcome **directly**, not through the metabolite | `$stage2_only` |
| $\mathcal G_4$ | $= 0$ | $= 0$ | not involved in the pathway | role `"none"` |

## Method

For a given reference taxon, estimation and inference take three steps:

1. **ALR transform.** Each taxon is divided by the reference taxon and logged, which removes the compositional
   constraint.
2. **Penalized two-stage least squares.**
   - *Stage 1* (external data): a lasso (Tibshirani 1996) of the metabolite on the ALR microbiome gives
     $\hat{\boldsymbol\gamma}$.
   - *Stage 2* (target data): the metabolite is predicted as $\hat m_i = \mathbf{x}_i^\top\hat{\boldsymbol\gamma}$,
     and a lasso of the outcome on $\hat m_i$ and the ALR microbiome gives $\hat\theta$ and $\hat{\boldsymbol\beta}$.
3. **Regularized permutation test** of $H_0:\theta^*=0$, adjusting for the microbes selected in stage 2
   (Hemerik et al. 2021).

**Robust multi-reference procedure.** Steps 1–3 are repeated with each of the $p$ microbes as the reference.
The per-reference p-values are combined with the Cauchy combination test (Liu and Xie 2020). The correlation
statistics are summarized by their median. Each microbe's selection frequency across references is then used for
stability selection (Meinshausen and Bühlmann 2010) with threshold $\tau \in (0.5, 1)$:

$$
\pi^{(\gamma)}_k = \frac{1}{p-1}\sum_{j\neq k} \mathbb 1\big(\hat\gamma^{(j)}_k \neq 0\big), \qquad
\pi^{(\beta)}_k = \frac{1}{p-1}\sum_{j\neq k} \mathbb 1\big(\hat\beta^{(j)}_k \neq 0\big).
$$

A microbe is assigned to $\mathcal G_1,\dots,\mathcal G_4$ according to whether $\pi^{(\gamma)}_k > \tau$ and/or
$\pi^{(\beta)}_k > \tau$.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("lei-fang-stat/MMDP")
```

The package depends on `glmnet` (Friedman et al. 2010) and `glmmTMB` (Brooks et al. 2017). `phd` is optional and
is used only in the tests.

## Main functions

| Purpose | Function | Returns |
|---|---|---|
| p-value for one metabolite (one reference) | `test_metabolite()` | p-value, method, $\hat\theta$, selected taxa |
| p-values for many metabolites / all references | `test_metabolites()`, `reference_sensitivity()` | data frame of p-values |
| Stage 1 estimation across references | `stage1_estimation()` | reference × taxon matrix of $\hat{\boldsymbol\gamma}$ |
| Stage 2 estimation across references | `stage2_estimation()` | reference × (correlation, taxa) matrix of $\hat{\boldsymbol\beta}$ |
| Stability selection of microbial roles | `stability_selection()` | microbes in $\mathcal G_2$ (stage 1 only), $\mathcal G_1$ (both), $\mathcal G_3$ (stage 2 only) |
| Helpers | `alr_transform()`, `prepare_microbiome()`, `flhd_test()`, `perm_cor_test()`, `simulate_cohorts()` | |

## Quick start (simulated data)

```r
library(predmetab)

sim <- simulate_cohorts(n_ref = 100, n_target = 80, n_taxa = 20, n_metab = 1)
m   <- sim$metabolites[, 1]

# 1. Test H0: theta* = 0 with a single reference taxon
prep <- prepare_microbiome(sim$ref_counts, sim$target_counts, reference = 20)
test_metabolite(m, prep$X_ref, prep$X_target, sim$y, nperm = 2000)

# 2. Robust multi-reference p-value: one test per reference, then Cauchy combination
pv <- reference_sensitivity(sim$ref_counts, sim$target_counts, sim$metabolites, sim$y,
                            nperm = 2000)
cauchy_combine <- function(p) {
  p <- pmin(pmax(p[!is.na(p)], 1e-15), 1 - 1e-15)
  0.5 - atan(mean(tan((0.5 - p) * pi))) / pi
}
cauchy_combine(pv$p_value)

# 3. Stage 1 and stage 2 estimates for every reference taxon
gamma_all <- stage1_estimation(m, sim$ref_counts, sim$target_counts)
theta_all <- stage2_estimation(m, sim$ref_counts, sim$target_counts, sim$y)
median(theta_all[, 1], na.rm = TRUE)   # median correlation across references

# 4. Stability selection of microbial roles (tau = 0.8)
roles <- stability_selection(gamma_all, theta_all, threshold = 0.8)
roles
#> Stability selection, threshold 0.8
#>   Stage 1 only (2): ...      <- G2: act through the metabolite
#>   Both stages  (2): ...      <- G1
#>   Stage 2 only (10): ...     <- G3: direct effect on the outcome
```

## Application: iHMP and FRAN (inflammatory bowel disease)

In the paper, the external dataset is iHMP (Lloyd-Price et al. 2019; 104 subjects with paired microbiome and
metabolomics). The target dataset is FRAN (Franzosa et al. 2019; 153 subjects with fecal calprotectin). Both use
86 shared microbial orders and 76 microbially related metabolites, from the curated collection of Muller et al.
(2022)(<https://github.com/borenstein-lab/microbiome-metabolome-curated-data/wiki/Data-overview#references>). A pseudo-abundance of $10^{-8}$ was added
before the ALR transform.

```r
IHMP.mtb   <- read.table("IHMP-mtb-main.txt",  sep = "\t")   # external metabolites
IHMP.micro <- read.table("IHMP.nozero.86.txt", sep = "\t")   # external microbiome
FRAN.micro <- read.table("FRAN.nozero.86.txt", sep = "\t")   # target microbiome
FRAN.meta  <- read.table("FRAN-Y.txt",         sep = "\t")   # target outcome

Y <- FRAN.meta[complete.cases(FRAN.meta), ]
colnames(Y) <- c("subject", "Fecal")
y <- as.numeric(scale(log(as.numeric(Y$Fecal)), scale = FALSE))   # centred log FCP
FRAN.micro.y <- FRAN.micro[match(Y$subject, rownames(FRAN.micro)), ]

# p-values for all metabolites, every microbe used as the reference
pv <- reference_sensitivity(IHMP.micro, FRAN.micro.y, IHMP.mtb, y,
                            references = 1:86, mc.cores = 20, nperm = 2e4)
p_combined <- tapply(pv$p_value, pv$metabolite, cauchy_combine)

# microbial roles for one metabolite
i <- which(colnames(IHMP.mtb) == "HMDB0000152")          # gentisic acid
gamma_all <- stage1_estimation(IHMP.mtb[, i], IHMP.micro, FRAN.micro)
theta_all <- stage2_estimation(IHMP.mtb[, i], IHMP.micro, FRAN.micro, y,
                               target_ids = Y$subject, m_name = colnames(IHMP.mtb)[i])
stability_selection(gamma_all, theta_all, threshold = 0.8)
```

## Citation

If you use `MMDP`, please cite:

> Fang L, Ye C, Wang Y. Integrating complementary datasets for constructing microbiome–metabolome–disease
> pathways. *Bioinformatics* (accepted).

## References
Lloyd-Price J, Arze C, Ananthakrishnan AN *et al.* Multi-omics of the gut microbial ecosystem in inflammatory bowel
diseases. *Nature* 2019;**569**:655–62.
Franzosa EA, Sirota-Madi A, Avila-Pacheco J *et al.* Gut microbiome structure and metabolic activity in
inflammatory bowel disease. *Nat Microbiol* 2019;**4**:293–305.
Muller E, Algavi YM, Borenstein E. The gut microbiome-metabolome dataset collection: a curated resource for
integrative meta-analysis. *NPJ Biofilms Microbiomes* 2022;**8**:79.

## License

GPL (>= 2)

## Contact

Lei Fang (fangl17@miamioh.edu) 
