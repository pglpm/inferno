data('iris')

devtools::install()
library('modelfreeinference')

buildmetadata(data = iris, file = 'iris_metadata.csv')
