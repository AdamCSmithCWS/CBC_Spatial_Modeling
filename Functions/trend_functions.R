
# helper functions --------------------------------------------------------

p_lt <- function(x,th){
  length(which(x < th))/length(x)
}


p_neg <- function(x){
  length(which(x < 0))/length(x)
}

texp <- function(x,ny = 2019-1974){
  (x^(1/ny)-1)*100
}


chng <- function(x){
  (x-1)*100
}

prob_dec <- function(ch,thresh){
  
  length(which(ch < thresh))/length(ch)
}
