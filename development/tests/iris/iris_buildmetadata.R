data('iris')

devtools::install()
library('predict')

buildmetadata(data = iris, file = '_test_iris_metadata.csv')
