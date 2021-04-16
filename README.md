# Extended functionalities for `rpart` [<img src="man/figures/rpart.png" alt="Rpart logo" style="float:right;height:232.25px" align="right" height="232.25">](https://cran.r-project.org/web/packages/rpart/index.html)

This modified version of the well-known [`rpart`](https://cran.r-project.org/web/packages/rpart/index.html) package extends the functionalities of the latter so that the user can:

- Modify the *mtry* parameter (number of predictors consider at each split).
- Use new split functions which allow for dealing with predictands that are non-normally distributed. For speed purposes, this part of the code has been written in the C language.

This allows for building non-standard **random forests** through the [RandomForest2](https://github.com/MNLR/RandomForest2) package. 

**Note**: This package is still being developed so new functionalities may be added soon.
