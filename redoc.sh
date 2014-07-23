#!/bin/sh

R --no-save <<EOF
library(devtools)
document()
EOF
