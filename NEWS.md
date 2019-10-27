### Changes in the cosa 2.1.0
 - Functions can accomodate polynomial functional forms for the score variable w/ or w/o interaction with treatment
 - Improvements to bound constrained optimization
 
### Changes in the cosa 2.0.0
 - Functions can accomodate quadratic functional forms for the score variable
 - Users can override default degrees of freedom
 - A few minor changes to naming conventions

### Changes in the cosa 1.2.2
 - Added `vectorize.cosa()` function, which makes it convenient to run and compare multiple designs
 - Removed MMA and COBYLA optimizers from default options
 
### Changes in the cosa 1.2.1
 - Removed experimental wrapper functions `cosa()`, `mdes()` and `power()`
 - Summaries are printed with the function call
 - Users can specify benchmark values in plots
 - Changed the package title

### Changes in the cosa 1.2.0
 - Added jacobians
 - Option for maximizing power when primary constraint is placed on the total cost
 - Bug fixes for `cosa.bird2r1()` function 

### Changes in the cosa 1.1.0
 - Added experimental wrapper functions `cosa()`, `mdes()` and `power()`
 - Added designs with top-level stratification (fixed blocks) 
 - Bug fixes for `bcrd4r2` and `bcrd4r3` designs
 - Added S3 `plot()` methods
 - Improved validity checks for arguments
 - Added citation
