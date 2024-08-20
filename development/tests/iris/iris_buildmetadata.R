data('iris')

devtools::install()
library('inferno')

buildmetadata(data = iris, file = '_test_iris_metadata.csv')
