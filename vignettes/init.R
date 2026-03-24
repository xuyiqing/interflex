
# List of packages
packages <- c("Rcpp", "mgcv", "sandwich", "pcse", "foreach",
              "doParallel", "lfe", "lmtest", "Lmoments","ggplot2",
              "plotrix", "grid", "gridExtra", "ggplotify", "ggpubr", "labelled",
              "RColorBrewer", "grDevices", "gtable", "MASS", "mvtnorm",
              "pROC", "ModelMetrics", "foreign","patchwork", "rmarkdown", "DT", "interflex")

# Install and load each package
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

