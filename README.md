# Extended functionalities for `rpart` [<img src="man/figures/rpart.png" alt="Rpart logo" style="float:right;height:232.25px" align="right" height="232.25">](https://cran.r-project.org/web/packages/rpart/index.html)

This modified version of the well-known [`rpart`](https://cran.r-project.org/web/packages/rpart/index.html) package, which has been forked from [here](https://github.com/bethatkinson/rpart), extends the functionalities of `rpart` so that the user can:

- Modify the *mtry* parameter (number of predictors considered at each split).
- Use new split functions which allow for dealing with predictands that are non-normally distributed. For speed purposes, this part of the code has been written in the C language.

This allows for building non-standard **random forests** through the [RandomForestDist](https://github.com/MNLR/RandomForestDist) package. 

**Note**: This package is still being developed and new functionalities may be added soon.

To install this [modified version of rpart](https://github.com/MNLR/rpart), simply type the following:

```
devtools::install_github("MNLR/rpart")
```

In case `devtools` is not already available, it can be installed from CRAN using the command `install.packages("devtools")`.



