# R code that I use to generate examples for teaching

This package is based on examples that I use to generate notes for various courses.  It is evolving into a holding tank for implementations of various population- and quantitative- genetic concepts.  It will likely never be a "real" R package in CRAN.

# Getting the code

```
git clone https://github.com/molpopgen/Rcode4teaching
```

# Installation

```
R CMD INSTALL Rcode4teaching
```

# Building the documentation

By "documentation", I mean those often-terse pdfs that come with R packages:

```
R CMD Rd2pdf Rcode4teaching
```
