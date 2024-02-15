# msDatabase <img src="man/figures/msDatabase.png" align="right" alt="msDatabase logo" style="height: 140px;"></a>

Tool for manipulating Open source mass spectrometry databases

## How to install?

``` r
install.packages("devtools")

devtools::install_github("Zhujiamin97/msDatabase", build_vignettes = TRUE)
```

## Run example
## [MASSBANK_EU](https://massbank.eu/MassBank/Search)
``` r
library(msDatabase)
result <- search.massbank_eu(compound_name = "Dihydrotestosterone",
                             formula = NULL,
                             InChIKey = NULL)
print(result)

# show spec
plot_spec(search.result = result,
          row_num = 2)

# match ms2
result <- search.massbank_eu(compound_name = "Dihydrotestosterone",
                             ms2_matrix = ms2_matrix_demo())
# get dp score
print(result$dp.score)
```
## [HMDB](https://hmdb.ca/)
``` r
library(msDatabase)
result <- search.hmdb(HMDB_ID = "HMDB0000001")
print(result)

# show spec
plot_spec(search.result = result,
          row_num = 2)
```
## [T3DB](http://www.t3db.ca/)
``` r
library(msDatabase)
result <- search.t3db(T3DB_ID = "T3D3216")
print(result)

# show spec
plot_spec(search.result = result,
          row_num = 2)
```

## [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
``` r
library(msDatabase)
result <- search.pubchem(cids = c(1,2,3))
print(result)
```
