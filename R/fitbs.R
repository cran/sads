fitbs <- function(x, trunc, ...){
    dots <-list(...)
    if (any(x <= 0)) stop ("All x must be positive")
    s <- length(x)
    n <- sum(x)
    if (!missing(x)){
        if (!missing(trunc)){
            if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
            else{
                LL <- function(N,S) -sum(dtrunc("bs", x = x, coef = list(N = N, S = S), trunc = trunc, log = TRUE))
            }
        }
        if (missing(trunc)){
            LL <- function(N,S) -sum(dbs(x = x, N = N, S = S, log = TRUE))
        }
        result <- do.call("mle2", c(list(minuslogl=LL, data = list(x = x), fixed=list(N=n, S=s), eval.only=TRUE), dots))
        new("fitsad", result, sad = "bs", distr = distr.depr, trunc = ifelse(missing(trunc), NaN, trunc))
    }
}
