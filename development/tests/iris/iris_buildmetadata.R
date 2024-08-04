data('iris')

devtools::install()
library('modelfreeinference')

buildmetadata(data = iris, file = '_test_iris_metadata.csv')
