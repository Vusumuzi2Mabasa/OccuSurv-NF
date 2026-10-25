# Install runtime dependencies into the library selected by R_LIBS_USER.
options(repos=c(CRAN='https://cloud.r-project.org'))
install.packages(c('survival','ggplot2','jsonlite','knitr','rmarkdown'))
