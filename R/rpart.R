#
#  The recursive partitioning function, for R
#
rpart <-
    function(formula, data, weights, subset, na.action = na.rpart, method,
             model = FALSE, x = FALSE, y = TRUE, parms, control, cost, 
             mtry = NA,
             ...)
{
    Call <- match.call()
    
    method.names <- c("anova", "poisson", "class", "exp", 
                      "gammaLLMME", "gammaLLmean", 
                      "bernoulliGammaLLMME",
                      "gammaDeviation", 
                      "gammaLLBC3",
                      "bernoulliLL",
                      "binaryCrossEntropyGammaDeviation",
                      "binaryCrossEntropy",
                      "binaryCrossEntropyMultivar",
                      "binaryCrossEntropyMultivarCorPenalization",
                      "binaryMultiEntropy",
                      "binaryMultiEntropyCond", 
                      "binaryMargEntropyCond", 
                      "multiBinaryGammaEntropy",
                      "MSEgammaDeviance",
                      "MSEbinaryEntropyGammaDeviance",
                      "binaryDoubleEntropyGammaDeviance", 
                      "multiMSE")
    


    # if (!missing(method) && method == "bernoulliLL" && 
    #     !all.equal(unname(sort(unique(y))), c(0,1)))
    #   stop("Either the variable has NAs or is not Bernoulli distributed")

    
    if (is.data.frame(model)) {
        m <- model
        model <- FALSE
    } else {
        indx <- match(c("formula", "data", "weights", "subset"),
                      names(Call), nomatch = 0L)
        if (indx[1] == 0L) stop("a 'formula' argument is required")
        temp <- Call[c(1L, indx)]      # only keep the arguments we wanted
        temp$na.action <- na.action    # This one has a default
        temp[[1L]] <- quote(stats::model.frame) # change the function called
        m <- eval.parent(temp)
    }

    Terms <- attr(m, "terms")
    if (any(attr(Terms, "order") > 1L))
	stop("Trees cannot handle interaction terms")

    wt <- model.weights(m)
    if (any(wt < 0)) stop("negative weights not allowed")
    if (!length(wt)) wt <- rep(1, nrow(m))
    offset <- model.offset(m)
    X <- rpart.matrix(m)
    
    
    nobs <- nrow(X)
    nvar <- ncol(X)
    
    Y <- model.response(m)
    
    if (method == "binaryMultiEntropy"){
      if ( !identical(unique(Y), c(0,1)) && !identical(unique(Y), c(1,0)) ){
        stop("Data is not binary") 
      }
    }
    
    if (method %in% c("binaryCrossEntropyMultivar",
                      "binaryMultiEntropy",
                      "binaryMultiEntropyCond",
                      "binaryMargEntropyCond",
                      "multiBinaryGammaEntropy" ,
                      "MSEgammaDeviance",
                      "MSEbinaryEntropyGammaDeviance",
                      "binaryDoubleEntropyGammaDeviance",
                      "multiMSE")){
        if (is.null(ncol(Y)) || ncol(Y) == 1) stop("For multivariable methods y must be a 2D matrix")
      parms <- ncol(Y)
    } else if ( method == "binaryCrossEntropyMultivarCorPenalization" ){
      if (is.null(ncol(Y)) || ncol(Y) == 1) stop("For method 'binaryCrossEntropyMultivarCorPenalization' y must be a 2D matrix")
      if ( !identical(unique(Y), c(0,1)) && !identical(unique(Y), c(1,0)) ){
        stop("Data is not binary") 
      }
      if (length(parms) != 3){
        warning("3 parameters were not provided for 'binaryCrossEntropyMultivarCorPenalization',
                 Using default cor = 'pearson'; alpha = 0.5; beta = 0.5")
        alpha <- 0.5
        beta <- 0.5
        cor.type <- "pearson"
      } else {
        cor.type <- match.arg(arg = parms[1], choices = c("spearman",
                                                          "pearson", 
                                                          "kendall"))
        alpha <- as.numeric(parms[2])
        beta <- as.numeric(parms[3])
        if (!is.numeric(alpha) || !is.numeric(beta) ||
            alpha < 0 || beta < 0) stop("provide valid numeric alpha, beta")
      }
      
      aux.cor_ <- cor(Y, method = cor.type)
      
      cor_ <- rep(x = 0, times = ((ncol(Y)-1)*ncol(Y))/2)
      
      icor_ <- 1
      for (ic in 1:ncol(aux.cor_)){
        
        ir <- 1
        while (ir <= nrow(aux.cor_)){
          if (ir > ic){
            cor_[icor_] <- aux.cor_[ir, ic]
            icor_ <- icor_ + 1
          }
          ir <- ir + 1
        }
      }
      
      parms <- c(ncol(Y), 
                 alpha, 
                 beta,
                 cor_)
    }
    
    #print("X = ")
    #print(X)
    #print("Y = ")
    #print(Y)
    
    if (is.null(mtry)){
      if (is.character(Y) || is.factor(Y)) mtry <- floor(sqrt(nvar))
      else mtry <- max(floor(nvar/3), 1)
    } else if (is.na(mtry)) mtry <- nvar
      else if ( x%%1 != 0 || mtry < 1 || mtry > nvar) stop("Provide a valid mtry")
    
    if (missing(method) || is.null(method)) { #ensures compatibility with wrapping
    	method <- if (is.factor(Y) || is.character(Y)) "class"
            else if (inherits(Y, "Surv")) "exp"
    	else if (is.matrix(Y)) "poisson"
    	else "anova"
    }

    if (is.list(method)) {
        ## User-written split methods
      mlist <- method
      method <- "user"

        ## Set up C callback.  Assign the result to a variable to avoid
        ## garbage collection
      init <- if (missing(parms) || is.null(parms)) mlist$init(Y, offset, wt = wt) else mlist$init(Y, offset, parms, wt)
      keep <- rpartcallback(mlist, nobs, init)
        
      method.int <- 4L             # the fourth entry in func_table.h
      #	numresp <- init$numresp
      #	numy <- init$numy
      parms <- init$parms
    } else {
      method.int <- pmatch(method, method.names)
      if (is.na(method.int)) stop("Invalid method")
      method <- method.names[method.int]
      if (method.int == 4L) method.int <- 2L

        ## If this function is being retrieved from the rpart package, then
        ##   preferentially "get" the init function from there.  But don't
        ##   lock in the rpart package otherwise, so that we can still do
        ##   standalone debugging.
	init <- if (missing(parms) || is.null(parms)) get(paste("rpart", method, sep = "."),
	                                envir = environment())(Y, offset, , wt)
        else
            get(paste("rpart", method, sep = "."),
                envir = environment())(Y, offset, parms, wt)
        ## avoid saving environment on fitted objects
        ns <- asNamespace("rpart")
        if (!is.null(init$print)) environment(init$print) <- ns
        if (!is.null(init$summary)) environment(init$summary) <- ns
        if (!is.null(init$text)) environment(init$text) <- ns
    }

    Y <- init$y

    xlevels <- .getXlevels(Terms, m)
    cats <- rep(0L, ncol(X))
    if (!is.null(xlevels))
	cats[match(names(xlevels), colnames(X))] <-
            unlist(lapply(xlevels, length))

    ## We want to pass any ... args to rpart.control, but not pass things
    ##  like "dats = mydata" where someone just made a typo.  The use of ...
    ##  is simply to allow things like "cp = 0.05" with easier typing
    extraArgs <- list(...)
    if (length(extraArgs)) {
	controlargs <- names(formals(rpart.control)) # legal arg names
	indx <- match(names(extraArgs), controlargs, nomatch = 0L)
	if (any(indx == 0L))
            stop(gettextf("Argument %s not matched",
                          names(extraArgs)[indx == 0L]),
                 domain = NA)
    }

    controls <- rpart.control(...)
    if (!missing(control) && !is.null(control)) controls[names(control)] <- control

    xval <- controls$xval
    if (is.null(xval) || (length(xval) == 1L && xval == 0L) || method=="user") {
	xgroups <- 0L
	xval <- 0L
    } else if (length(xval) == 1L) {
        ## make random groups
        xgroups <- sample(rep(1L:xval, length = nobs), nobs, replace = FALSE)
    } else if (length(xval) == nobs) {
	xgroups <- xval
	xval <- length(unique(xgroups))
    } else {
        ## Check to see if observations were removed due to missing
	if (!is.null(attr(m, "na.action"))) {
            ## if na.rpart was used, then na.action will be a vector
	    temp <- as.integer(attr(m, "na.action"))
	    xval <- xval[-temp]
	    if (length(xval) == nobs) {
		xgroups <- xval
		xval <- length(unique(xgroups))
            } else stop("Wrong length for 'xval'")
        } else stop("Wrong length for 'xval'")
    }

    ##
    ## Incorporate costs
    ##
    if (missing(cost) || is.null(control)) cost <- rep(1, nvar)
    else {
	if (length(cost) != nvar)
            stop("Cost vector is the wrong length")
	if (any(cost <= 0)) stop("Cost vector must be positive")
    }

    ##
    ## Have C code consider ordered categories as continuous
    ##  A right-hand side variable that is a matrix forms a special case
    ## for the code.
    ##
    tfun <- function(x)
	if (is.matrix(x)) rep(is.ordered(x), ncol(x)) else is.ordered(x)
    labs <- sub("^`(.*)`$", "\\1", attr(Terms, "term.labels")) # beware backticks
    isord <- unlist(lapply(m[labs], tfun))

    storage.mode(X) <- "double"
    storage.mode(wt) <- "double"
    temp <- as.double(unlist(init$parms))
    if (!length(temp)) temp <- 0    # if parms is NULL pass a dummy
    rpfit <- .Call(C_rpart,
                   ncat = as.integer(cats * !isord),
                   method = as.integer(method.int),
                   as.double(unlist(controls)),
                   temp,
                   as.integer(xval),
                   as.integer(xgroups),
                   as.double(t(init$y)),
                   X,
                   wt,
                   as.integer(init$numy),
                   as.double(cost), 
                   as.integer(mtry))

    nsplit <- nrow(rpfit$isplit) # total number of splits, primary and surrogate
    ## total number of categorical splits
    ncat <- if (!is.null(rpfit$csplit)) nrow(rpfit$csplit) else 0L
#    nodes <- nrow(rpfit$inode)
    if (nsplit == 0L) xval <- 0L # No xvals were done if no splits were found

    numcp <- ncol(rpfit$cptable)
    temp <- if (nrow(rpfit$cptable) == 3L) c("CP", "nsplit", "rel error")
            else c("CP", "nsplit", "rel error", "xerror", "xstd")
    dimnames(rpfit$cptable) <- list(temp, 1L:numcp)

    tname <- c("<leaf>", colnames(X))
    splits <- matrix(c(rpfit$isplit[, 2:3], rpfit$dsplit), ncol = 5L,
                     dimnames = list(tname[rpfit$isplit[, 1L] + 1L],
                     c("count", "ncat", "improve", "index", "adj")))
    index <- rpfit$inode[, 2L]  # points to the first split for each node

    ## Now, make ordered factors look like factors again (a printout choice)
    nadd <- sum(isord[rpfit$isplit[, 1L]])
    if (nadd > 0L) { # number of splits at an ordered factor.
	newc <- matrix(0L, nadd, max(cats))
	cvar <- rpfit$isplit[, 1L]
	indx <- isord[cvar]             # vector of TRUE/FALSE
	cdir <- splits[indx, 2L]        # which direction splits went
	ccut <- floor(splits[indx, 4L]) # cut point
	splits[indx, 2L] <- cats[cvar[indx]] # Now, # of categories instead
	splits[indx, 4L] <- ncat + 1L:nadd # rows to contain the splits

        ## Next 4 lines can be done without a loop, but become indecipherable
	for (i in 1L:nadd) {
	    newc[i, 1L:(cats[(cvar[indx])[i]])] <- -as.integer(cdir[i])
	    newc[i, 1L:ccut[i]] <- as.integer(cdir[i])
        }
	catmat <- if (ncat == 0L) newc
        else {
            ## newc may have more cols than existing categorical splits
            ## the documentation says that levels which do no exist are '2'
            ## and we add 2 later, so record as 0 here.
            cs <- rpfit$csplit
            ncs <- ncol(cs); ncc <- ncol(newc)
            if (ncs < ncc) cs <- cbind(cs, matrix(0L, nrow(cs), ncc - ncs))
            rbind(cs, newc)
        }
	ncat <- ncat + nadd
    } else catmat <- rpfit$csplit

    ## NB: package adabag depends on 'var' being a factor.
    if (nsplit == 0L) {                    # tree with no splits
	frame <- data.frame(row.names = 1L,
			    var = "<leaf>",
			    n = rpfit$inode[, 5L],
			    wt = rpfit$dnode[, 3L],
			    dev = rpfit$dnode[, 1L],
			    yval = rpfit$dnode[, 4L],
			    complexity = rpfit$dnode[, 2L],
			    ncompete = 0L,
			    nsurrogate = 0L)
    } else {
	temp <- ifelse(index == 0L, 1L, index)
	svar <- ifelse(index == 0L, 0L, rpfit$isplit[temp, 1L]) # var number
	frame <- data.frame(row.names = rpfit$inode[, 1L],
                            ## maybe better to specify tname as the level?
			    var = tname[svar + 1L],
			    n = rpfit$inode[, 5L],
			    wt = rpfit$dnode[, 3L],
			    dev = rpfit$dnode[, 1L],
			    yval = rpfit$dnode[, 4L],
			    complexity = rpfit$dnode[, 2L],
			    ncompete = pmax(0L, rpfit$inode[, 3L] - 1L),
			    nsurrogate = rpfit$inode[, 4L])
    }
    if (method.int == 3L) {
        ## Create the class probability vector from the class counts, and
        ##   add it to the results
        ## Also scale the P(T) result
        ## The "pmax" 3 lines down is for the case of a factor y which has
        ##   no one at all in one of its classes.  Both the prior and the
        ##   count will be zero, which led to a 0/0.
        numclass <- init$numresp - 2L
        nodeprob <- rpfit$dnode[, numclass + 5L] / sum(wt) # see ginidev.c
        temp <- pmax(1L, init$counts)   # overall class freq in data
        temp <- rpfit$dnode[, 4L + (1L:numclass)] %*% diag(init$parms$prior/temp)
        yprob <- temp /rowSums(temp)    # necessary with altered priors
        yval2 <- matrix(rpfit$dnode[, 4L + (0L:numclass)], ncol = numclass + 1L)
	frame$yval2 <- cbind(yval2, yprob, nodeprob)
    } else if (init$numresp > 1L)
        frame$yval2 <- rpfit$dnode[, -(1L:3L), drop = FALSE]

    if (is.null(init$summary))
        stop("Initialization routine is missing the 'summary' function")
    functions <- if (is.null(init$print)) list(summary = init$summary)
                 else list(summary = init$summary, print = init$print)
    if (!is.null(init$text)) functions <- c(functions, list(text = init$text))
    if (method == "user") functions <- c(functions, mlist)

    where <- rpfit$which
    names(where) <- row.names(m)

    ans <- list(frame = frame,
             where = where,
             call = Call, terms = Terms,
             cptable = t(rpfit$cptable),
             method = method,
             parms = init$parms,
             control = controls,
             functions = functions,
             numresp = init$numresp)
    if (nsplit) ans$splits = splits
    if (ncat > 0L) ans$csplit <- catmat + 2L
    if (nsplit) ans$variable.importance <- importance(ans)
    if (model) {
	ans$model <- m
	if (missing(y) || is.null(y)) y <- FALSE
    }
    if (y) ans$y <- Y
    if (x) {
	ans$x <- X
	ans$wt <- wt
    }
    ans$ordered <- isord
    if (!is.null(attr(m, "na.action"))) ans$na.action <- attr(m, "na.action")
    if (!is.null(xlevels)) attr(ans, "xlevels") <- xlevels
    if (method == "class") attr(ans, "ylevels") <- init$ylevels
    class(ans) <- "rpart"
    ans
}
