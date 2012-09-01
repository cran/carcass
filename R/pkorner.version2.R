pkorner.version2 <-
function(s, f, d, n, k=0.25){
# s = probability that a dead bat is still there after 24 hours.
# f = searchers efficiency, probability that a dead bat is detected by a searcher
# d = (average) number of days between two searches
# n = number of searches (n * d = study period)
# k = factor by which serachers efficiency decreases per search, default = 0.25 as assumed by Manuela Huso                                  
A<-s*(1-s^d)/(1-s)
summepfound<-A*f
for(x in 2:n){
      summep<-1+k*s^d*(1-f)
      for(j in 1:(x-1)) summep<-summep+k^(x-j) * s^((x-j)*d) * productsef(f, k, n=x, j)
      summepfound<-summepfound+A*f*summep
      }
p<-summepfound/(d*n) 
p
}

