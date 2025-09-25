pkgs <- readLines("R-requirements.txt")
# cran <- readLines("R-requirements-cran.txt")
# bioc <- readLines("R-requirements-bioc.txt")

if (!require("BiocManager")) install.packages("BiocManager")

#Checking if the package belongs to CRAN or to Bioconductor and installing them accordingly.

for(lib in pkgs){
        if(!lib %in% installed.packages()){
            if(lib %in% available.packages()[,1]){
             install.packages(lib,dependencies=TRUE)
            }else {(BiocManager::install(lib))
            }}
        }

#Loading the libraries
sapply(pkgs,require,character=TRUE)