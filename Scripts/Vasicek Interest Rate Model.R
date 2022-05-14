## Simulate Sample Paths ##

## define model parameters
r0 <- 0.03                                                                #initial yield
theta <- 0.10                                                             #long term yield
k <- 0.3                                                                  #speed of reversion to long term yield
beta <- 0.03                                                              #instantaneous volatility

## simulate short rate paths
n <- 10    # MC simulation trials
T <- 10    # total time
m <- 200   # subintervals
dt <- T/m  # difference in time each subinterval

VasicekZCBprice_simulated <-
  function(num_trials, num_subintervals, total_time, r0, theta, k, beta)      #define function VasicekZCBprice_simulated
  {
    dt <- total_time/num_subintervals                                         #calculate dt to be the total time divided by number of intervals
    r <- matrix(0,num_subintervals+1,num_trials) # matrix to hold short rate paths           #create matrix with m+1 rows and n columns with each cell containing 0
    r[1,] <- r0                                                               #set all n values in row 1 to r0
    
    for(j in 1:num_trials){                                                   #for each sample path
      for(i in 2:(num_subintervals+1)){                                       #for each row starting at row 2
        Wt <- rnorm(1,0,1)                                                    #select 1 random number using standard normal distrubtion to use as Brownian motion component
        dr <- k*(theta-r[i-1,j])*dt + beta*sqrt(dt)*Wt                        #change in rate achieved using Vasicek model
        r[i,j] <- r[i-1,j] + dr                                               #new rate = previous rate plus change in rate; save to i'th row and j'th column of returns matrix
      }
    }
    return(r)                                                                 #return the create matrix r
  }

returns = VasicekZCBprice_simulated(n,m,T,r0,theta,k,beta)

## plot paths
t <- seq(0, T, dt)                                                        #create array of time steps starting at 0 and ending at T, counting by dt
rT.expected <- theta + (r0-theta)*exp(-k*t)                               #create array of points expected values at times t
rT.stdev <- sqrt( beta^2/(2*k)*(1-exp(-2*k*t)))                           #create array containing instantaneous standard deviations
matplot(t, returns[,1:10], type="l", lty=1, main="Short Rate Paths", ylab="rt") #plot each column of returns matrix 
abline(h=theta, col="red", lty=2)                                         #add red horizontal line at theta value 
abline(h=0)                                                               #add horizontal line at 0
lines(t, rT.expected, lty=2)                                              #add line of expected values to plot
lines(t, rT.expected + 2*rT.stdev, lty=2)                                 #add line of expected values + 2 standard deviations
lines(t, rT.expected - 2*rT.stdev, lty=2)                                 #add line of expected values - 2 standard deviations
points(0,r0)                                                              #add a point at (0,r0)

## function to find ZCB price using Vasicek model
VasicekZCBprice_exact <-                                                  #define function named VasicekZCBprice_exact that takes arguments r0, k, theta, beta, and T 
  function(r0, k, theta, beta, T){
    b.vas <- (1/k)*(1-exp(-T*k))                                          #within function calculate intermediate value and assign to variable 'b'
    a.vas <- (theta-beta^2/(2*k^2))*(T-b.vas)+(beta^2)/(4*k)*b.vas^2      #within function calculate intermediate value and assign to variable 'a'
    return(exp(-a.vas-b.vas*r0))                                          #return a value of e^(-a-b*r0) from VasicekZCBprice_exact function
  }

r0 <- seq(0.00, 0.20, 0.05)                                               #create a sequence counting from 0 to 0.2 by 0.05 and assign to r0
n <- length(r0)                                                           #set n to be the number of elements in the new r0 sequence
maturity <- seq(1, 10, 1)                                                 #create an array of times from 1 to 10 stepping by 1
m <- length(maturity)                                                     #set m to the number of values in maturity array
yield <- matrix(0, m, n)                                                  #create a m row x n column matrix with all elements having a value of 0
for(i in 1:n){                                                            #for each i from 1 to n
  for(T in 1:m){                                                          #for each T in sequence 1 to 10
    yield[T,i] <- -log(VasicekZCBprice_exact(r0[i], k, theta, beta, T))/T #set each value at row T, col i in yield matrix
  }
}

matplot(maturity, yield, type="l", col="black", lty=1, main="Yield Curve Shapes")  #plot yield curves in a line plot named 'Yield Curve Shapes' and using maturity values for x-axis
abline(h=theta, col="red", lty=2)                                         #add horizontal red line at theta

## simulate short rate paths
r0 <- 0.03     
n <- 10000  # MC simulation trials
T <- 1     # total time
m <- 200   # subintervals
dt <- T/m  # difference in time each subinterval

r = VasicekZCBprice_simulated(n,m,T,r0,theta,k,beta)                      #call VasicekZCBprice_simulated function

## calculate Monte Carlo bond price and compare to Exact Vasicek solution
ss <- colSums(r[2:(m+1),]*dt)  # integral estimate                        #create weighted average rate on each of the simulations
c <- exp(-ss)                                                             #calculate the bond price using each of the weighted average rates
estimate <- mean(c)                                                       #find the average of the prices found through simulations
stdErr <- sd(c)/sqrt(n)                                                   #calculate the stdErr of prices using standard deviation number of simulations
exact <- VasicekZCBprice_exact(r0, k, theta, beta, T)                     #calculate an exact price using VasicekZCBprice_exact function

                                                                          #display exact, estimated, and stdErr values (insert newlines '\n' for better readability)
cat('\n', 'Exact Vasicek Price:', round(exact,6), '\n',
    'MC Price:', round(estimate,6), '\n',
    'MC Standard Error:', round(stdErr,5))