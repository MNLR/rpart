# The `rpart` package <img src="man/figures/rpart.png" alt="Rpart logo" style="float:right;height:232.25px" align="right" height="232.25">


This is the source code for a modification of the [`rpart` package](https://cran.r-project.org/web/packages/rpart/index.html) — The recomended Classification and Regression Trees package in R — aimed to be used as the base tree package for building **random forests** by means of the package [RandomForest2](https://github.com/MNLR/RandomForest2).


## Overview

The `rpart` code builds classification or regression models of a very
general structure using a two stage procedure; the resulting models can be
represented as binary trees. The package implements many of the ideas found
in the CART (Classification and Regression Trees) book and programs of
Breiman, Friedman, Olshen and Stone.  Because CART is the trademarked name
of a particular software implementation of these ideas and `tree` was used
for the Splus routines of Clark and Pregibon, a different acronym -
Recursive PARTitioning or rpart - was chosen.
