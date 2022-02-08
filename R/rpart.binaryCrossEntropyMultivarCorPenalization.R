rpart.binaryCrossEntropyMultivarCorPenalization <- function(y, offset, parms, wt)
{
  if (!is.null(offset)) y <- y - offset
  list(y = y, parms = parms, numresp = ncol(y), numy = ncol(y),
       summary = function(yval, dev, wt, ylevel, digits ) {
         paste0("  mean=", formatg(yval, digits),
                ", err=" , formatg(dev, digits))
       },
       text = function(yval, dev, wt, ylevel, digits, n, use.n ) {
         if (use.n) paste0(formatg(yval, digits), "\nn=", n) else
           formatg(yval, digits)
       })
}
