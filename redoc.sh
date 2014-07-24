#!/bin/sh

R --no-save <<EOF
library(Rcpp)
library(devtools)
compileAttributes()
document()
EOF
