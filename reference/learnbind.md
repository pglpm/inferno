# Bind 3D arrays by first dimension

Used in 'util_checkpoints()' within 'learn()', and in various functions
in 'util_lprobs.R'.

## Usage

``` r
learnbind(x, y)
```

## Details

NB: the following variant is slower:

    function(x, y) {
        out <- c(aperm(x), aperm(y))
        dim(out) <- c(rev(dim(x)[-1]), dim(x)[1] + dim(y)[1])
        aperm(out)
    }
