# msDatabase <img src="man/figures/msDatabase.png" align="right" alt="msDatabase logo" style="height: 140px;"></a>

Tool for manipulating Open source mass spectrometry databases

## How to install?

``` r
install.packages("devtools")

devtools::install_github("Zhujiamin97/msDatabase", build_vignettes = TRUE)
```

## run example
``` r
library(msDatabase)
result <- search.massbank_eu(compound_name = "Dihydrotestosterone",
                             formula = NULL,
                             InChIKey = "QADHLRWLCPCEKT-LOVVWNRFSA-N")
print(result)

# show spec
plot_spec(search.result = result,
          row_num = 2)
```
