# Bound Constrained Optimal Sample Size Allocation (BCOSSA)

To install and load the package:
```{r}
install.packages("cosa")
library(cosa)
```

A limited version of BCOSSA is implemented in the shiny app (e.g., equal cost, fixed `p`), along with power and MDES functions. Live shiny app at: <br>
<https://cosa.shinyapps.io/index/> <br>

**cosa** implements bound constrained optimal sample size allocation (BCOSSA) framework described in Bulus & Dong (2019) for multilevel regression discontinuity designs (MRDDs) and multilevel randomized trials (MRTs) with continuous outcomes. It also includes functions to compute power rate or minimum detectable effect size. BCOSSA functions are designed to optimize proportion of treatment allocation (`p`) and sample size (generically, `n`) at one or more levels subject to budget, statistical power, or effect size constraints along with constraints on `n` and `p`. Constraints on `n` and `p` can be in the form of fixed values or bound constraints (aka box constraints).

Specifying `order = 0` or `rhots = 0` produces result equivalent to the corresponding random assignment design where also `p` can be optimized. `rhots = 0` means there is no relationship between the treatment [random] and the score variable, or `order = 0` means no score variable is included in the correctly specified model, therefore it is equavalent to random assignment designs.

Note: `n` and `p` should be omitted (or specified as `NULL`) for optimization. `p` can only be optimized (or bound constrained) in MRTs when treatment and control units have differing costs. When primary constraint is placed on a power rate or an effect size, providing marginal cost information will produce a cost-efficient sample allocation. When marginal cost information is not provided and two or more parameters are optimized, different starting values and algorithms may produce different results because the design is not uniquely identified. In such cases constraint `p = .50` is placed for MRTs. Comparing several algorithms and trying different starting values may faciliate decisions regarding sample sizes and `p` in such cases. One way is to set starting values at the expected values, or specify bounds such that they cover expected values. 


```{r}
# linear form interacting with the treatment
score.obj <- inspect.score(rnorm(1000), cutoff = 0,
                           order = 1, interaction = TRUE)

## single site (no blocks)
power.crd2(score.obj,
           es = .25, rho2 = .20, g2 = 0, r22 = 0,
           n1 = 50, n2 = 15)

## multiple sites (10 blocks)
## note that r22 > 0 due to explanatory power of indciator variables for sites
power.bcrd3f2(score.obj, 
              es = .25, rho2 = .20, g2 = 0, r22 = .30,
              n1 = 50, n2 = 15, n3 = 10)

## minimum required number of level 2 units per site
cosa.bcrd3f2(score.obj, 
             rho2 = .20, g2 = 0, r22 = .30,
             n1 = 50, n2 = NULL, n3 = 10)


# quadratic form interacting with the treatment
score.obj2 <- inspect.score(rnorm(1000), cutoff = 0,
                           order = 2, interaction = TRUE)

## minimum required number of level 2 units per site
cosa.bcrd3f2(score.obj2, 
             order = 2, interaction = TRUE, 
             rho2 = .20, g2 = 0, r22 = .30,
             n1 = 50, n2 = NULL, n3 = 10)
```

**Suggested citation**:

Bulus, M. (2021). Minimum Detectable Effect Size Computations for Cluster-Level Regression Discontinuity: Specifications Beyond Linear Functional Form. *Journal of Research on Educational Effectiveness*. Conditionally accepted.

Bulus, M., & Dong, N. (2019). Bound Constrained Optimization of Sample Sizes Subject to Monetary Restrictions in Planning of Multilevel Randomized Trials and Regression Discontinuity Studies. *The Journal of Experimental Education*. Advance online publication. <https://doi.org/10.1080/00220973.2019.1636197>

--o-- 
