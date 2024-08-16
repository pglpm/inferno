data(mtcars)

devtools::install()
library('predict')

buildmetadata(data = mtcars, file = '_test_mtcars_metadata.csv')
