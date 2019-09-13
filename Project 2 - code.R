library(ggplot2)

# Exercise a)
tau = 0.4
x = c(sample(c(0,1),1),rep(NA,249))
y = c(x[1] + rnorm(1,0,tau), rep(NA,249))
for (i in 2:250){
  if (runif(1)>0.1){
    x[i] = x[i-1]
  }
  else{x[i] = 1-x[i-1]}
  y[i] = x[i] + rnorm(1,0,tau)
}
df = data.frame(x = x, y = y)
ggplot(df, aes(x = 1:250, y = x)) + geom_line() + labs(x = "x", title = "Simulated values for x") + theme_classic()