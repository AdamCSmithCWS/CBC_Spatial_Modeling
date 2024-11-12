# data
library(mgcv)
ny <- 54
x <- 1:ny
fx <- 10 + 3*sin(x/10) + (-0.03*x)
y <- rpois(ny, fx)
x <- x + 1965
y <- log(y)
kts <- round(ny/4)
plot(x=x, y=exp(y), type="l")
lines(y=exp(fitted(gam(y ~ 1 + s(x, k=kts, bs="ps"), family="gaussian"))),
      x=x)

# gam trend function. p-spline with nyears/4 knots. gaussian on log Y
tgam <- function(y, x, t1, t2){
  require(mgcv)
  ny <- length(x) # full time series
  kts <- as.numeric(round(ny/4, 0)) # knot every four years
  df <- data.frame(y=y, x=x)
  # predict from p-spline with ny/4 knots, gaussian on log Y
  yhat <- exp(fitted(gam(y ~ 1 + s(x, k=kts, bs="ps"),
                         family="gaussian", data=df)))
  nt1 <- yhat[which(x==t1)] # pick start year
  nt2 <- yhat[which(x==t2)] # pick end year
  trd <- 100 * (((nt2 / nt1)^(1 / (t2 - t1))) - 1) # calc annl percent chg
  return(trd)
}

# test
t1 <- 1970
t2 <- 2019
nyrs <- t2 - t1 + 1
t <- tgam(y=y, x=x, t1=t1, t2=t2); t # annual percent change from gam endpoints
chg <- ((((t/100)+1)^nyrs)-1)*100; chg # accumulated percent change from gam endpoints
