# TreeSS
A model-free subdata selection method which uses CART-type trees to find a subdata which could be used for prediction purposes

## Development status

This package is currently being developed.

## Installation


The latest version of the package under development can be installed from GitHub:

```{r}
install.packages("devtools")
library(devtools)
remotes::install_github("agrakhi/TreeSS")
```

### Bug reports

Please submit any bugs or issues (or suggestions) using the [issues](https://github.com/agrakhi/TreeSS/issues) tab of the repo.

## Usage

The main functions users will use are `trees_subdata_withTwin` and `trees_subdata_withUni`. The former selects sample based on Twinning, while the latter does so based on uniform samples. Both use tree-based partition before sampling.

Check out the vignettes for more examples and details.

## License

This package is released in the public domain under the General Public License [GPL](https://www.gnu.org/licenses/gpl-3.0.en.html). 

## References

Singh, R., Stufken J.. “Model-free tree-based subdata selection.” doi: []().



