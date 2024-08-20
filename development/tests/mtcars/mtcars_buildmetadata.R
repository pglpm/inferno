data(mtcars)

devtools::install()
library('inferno')

buildmetadata(data = mtcars, file = '_test_mtcars_metadata.csv')
