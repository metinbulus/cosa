# Bound Constrained Optimal Sample Allocation (BCOSA)

To install and load the package:
```{r}
install.packages("cosa")
library(cosa)
```

A limited version of BCOSA is implemented in the shiny app (e.g., equal cost, fixed `p`), along with power and MDES functions. Live shiny app at: <br>
<https://cosa.shinyapps.io/index/> <br>

**cosa** implements bound constrained optimal sample size allocation (BCOSSA) framework described in Bulus & Dong (2019) for multilevel regression discontinuity designs (MRDDs) and multilevel randomized trials (MRTs) with continuous outcomes. It also includes functions to compute power rate or minimum detectable effect size. BCOSSA functions are designed to optimize proportion of treatment allocation (`p`) and sample size (generically, `n`) at one or more levels subject to budget, statistical power, or effect size constraints along with constraints on `n` and `p`. Constraints on `n` and `p` can be in the form of fixed values or bound constraints (aka box constraints).

Specifying `order = 0` or `rhots = 0` produces result equivalent to the corresponding random assignment design where also `p` can be optimized. `rhots = 0` means there is no relationship between the treatment [random] and the score variable, or `order = 0` means no score variable is included in the correctly specified model, therefore it is equavalent to random assignment designs.

Note: `n` and `p` should be omitted (or specified as `NULL`) for optimization. `p` can only be optimized (or bound constrained) in MRTs when treatment and control units have differing costs. When primary constraint is placed on a power rate or an effect size, providing marginal cost information will produce a cost-efficient sample allocation. When marginal cost information is not provided and two or more parameters are optimized, different starting values and algorithms may produce different results because the design is not uniquely identified. In such cases constraint `p = .50` is placed for MRTs. Comparing several algorithms and trying different starting values may faciliate decisions regarding sample sizes and `p` in such cases. One way is to set starting values at the expected values, or specify bounds such that they cover expected values. 


```{r}
score.obj <- inspect.score(rnorm(10000), cutoff = 0,
                           order = 2, interaction = TRUE)
                           
power.crd2(score.obj,
           es = .25, rho2 = .20, g2 = 0, r22 = 0,
           n1 = 50, n2 = 100)

# with 5 blocks df = n2 - 2*(n blocks) - order * (1 + int) - g2
# int  = 1 if interaction = TRUE
# n2: number of level 2 units across five blocks

power.crd2(score.obj, df = 100 - 2*5 - 2 * (1 + 1) - 1,
           es = .25, rho2 = .20, g2 = 0, r22 = .30,
           n1 = 50, n2 = 100)

# compare
# n2: number of level 2 units per block, n3: number of blocks

power.bcrd3f2(score.obj, 
              es = .25, rho2 = .20, g2 = 0, r22 = .30,
              n1 = 50, n2 = 20, n3 = 5)

# optimal combination of sample sizes for level 1 and level 2
# that produce power = .80 (given range restriction for level 1 sample size)

score.var <- rnorm(10000)

# second order polynomial order

cosa.bcrd3f2(score.var, cutoff = .50, treat.lower = TRUE,
             order = 2, interaction = FALSE,
             constrain = "power", power = .80,
             rho2 = .20, g2 = 0, r22 = .30,
             n1 = c(20, 60), n2 = NULL, n3 = 5)
             
# second order polynomial order interacting with treatment

cosa.bcrd3f2(score.var, cutoff = .50, treat.lower = TRUE,
             order = 3, interaction = TRUE, 
             constrain = "power", power = .80,
             rho2 = .20, g2 = 0, r22 = .30,
             n1 = c(20, 60), n2 = NULL, n3 = 5)
```

**Suggested citation**:

Bulus, M., & Dong, N. (2019). Bound Constrained Optimization of Sample Sizes Subject to Monetary Restrictions in Planning of Multilevel Randomized Trials and Regression Discontinuity Studies. *The Journal of Experimental Education*. Advance online publication. <https://doi.org/10.1080/00220973.2019.1636197>

--o-- 
