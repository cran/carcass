pkorner <-
function(s, f, d, n, k=0.25, search.efficiency.constant=TRUE){
# s = probability that a carcass remains 24 hours
# f = probability that a carcass is detected by a 
# searcher during a search given it persisted to the search
# d = (average) number of days between two searches
# n = number of searches (n * d = length of study period)
#--------------------------------------------------------------
if(search.efficiency.constant){
x <- (1-f)*s^d
A <- s*(1-s^d)/(1-s)
summep <- numeric(n)
for(k in 0:(n-1)) summep[k+1] <- (n-k)*x^k
p <- A*f*sum(summep)/(d*n)
}
if(!search.efficiency.constant){
p <- pkorner.version2(s=s, f=f, d=d, n=n, k=k)
}
p
}

