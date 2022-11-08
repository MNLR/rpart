rpart.multiMSE <- function(y, offset, parms, wt)
{
  if (!is.null(offset)) y <- y - offset
  list(y = y, parms = NULL, numresp = 2L, numy = 2L,
       summary = function(yval, dev, wt, ylevel, digits ) {
         paste0("  mean=", formatg(yval, digits),
                ", RMSE=" , formatg(dev, digits))
       },
       text = function(yval, dev, wt, ylevel, digits, n, use.n ) {
         if (use.n) paste0(formatg(yval, digits), "\nn=", n) else
           formatg(yval, digits)
       })
}
