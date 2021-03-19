# `rpart` for Random Forests with more Split Functions


This is the source code for a modification of the [`rpart` package](https://cran.r-project.org/web/packages/rpart/index.html) — The recomended Classification and Regression Trees package in R — aimed to be used as the base tree package for building **random forests** by means of the package [RandomForest2](https://github.com/MNLR/RandomForest2).

Essentially, this forked package has two main modifications:

- The internal code has been modified so that building trees now accepts the parameter mtry.
- New split functions have been added, mainly to deal with predictands that are gamma distributed. These are coded internally in c. The original package does contain an extension mechanism for R, but using split functions in R is extremelly slow. Note that this is work in progress.

## Overview of the original `rpart` Package [<img src="man/figures/rpart.png" alt="Rpart logo" style="float:right;height:232.25px" align="right" height="232.25">](https://cran.r-project.org/web/packages/rpart/index.html)

The `rpart` code builds classification or regression models of a very
general structure using a two stage procedure; the resulting models can be
represented as binary trees. The package implements many of the ideas found
in the CART (Classification and Regression Trees) book and programs of
Breiman, Friedman, Olshen and Stone.  Because CART is the trademarked name
of a particular software implementation of these ideas and `tree` was used
for the Splus routines of Clark and Pregibon, a different acronym -
Recursive PARTitioning or rpart - was chosen.
