data(mtcars)

devtools::install()
library('modelfreeinference')

buildmetadata(data = mtcars, file = '_test_mtcars_metadata.csv')
