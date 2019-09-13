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

ggplot(df, aes(x = 1:250, y = y)) + geom_point() + labs(x = "x", title = "Simulated values for y") + theme_classic()

# Exercise b)

likelihood <- function(parameters){
  p = parameters[1]
  tau = parameters[2]
  # First initiate normality constants and forward probabilities
  norm_const = c(1/(dnorm(y[1],0,tau)*0.5+dnorm(y[1],1,tau)*0.5), rep(NA,249))
  forward_prob = matrix(0,nrow=250,ncol = 2)
  forward_prob[1,1] = norm_const[1]*0.5*dnorm(y[1],0,tau)
  forward_prob[1,2] = norm_const[1]*0.5*dnorm(y[1],1,tau)
  transition = matrix(NA, nrow = 250, ncol = 4)
  
  # Iterate to find constants and probabilities
  for (t in 2:250){
    norm_const[t] = 1/(dnorm(y[t],0,tau)*p*forward_prob[t-1,1] + (1-p)*dnorm(y[t],0,tau)*forward_prob[t-1,2] + dnorm(y[t],1,tau)*(1-p)*forward_prob[t-1,1] + dnorm(y[t],1,tau)*p*forward_prob[t-1,2])
    transition[t,1] = p*dnorm(y[t],0,tau)*forward_prob[t-1,1]*norm_const[t]
    transition[t,2] = (1-p)*dnorm(y[t],1,tau)*forward_prob[t-1,1]*norm_const[t]
    transition[t,3] = (1-p)*dnorm(y[t],0,tau)*forward_prob[t-1,2]*norm_const[t]
    transition[t,4] = p*dnorm(y[t],1,tau)*forward_prob[t-1,2]*norm_const[t]
    forward_prob[t,1] = transition[t,1] + transition[t,3]
    forward_prob[t,2] = transition[t,2] + transition[t,4]
  }
  return(log(1/prod(norm_const)))
}

neg_likelihood <- function(parameters){
  return(-likelihood(parameters))
}

grid_points = 50
p_values = rep(seq(0, 1, length.out = grid_points),grid_points)
tau_values = rep(seq(0.1, 1, length.out = grid_points),each = grid_points)

test = matrix(c(p_values, tau_values), ncol = 2, byrow = F)
likelihood_values = matrix(apply(test,1,likelihood), ncol= grid_points)

x_axis = unique(p_values)
y_axis = unique(tau_values)

plot_data <- list(x_axis, y_axis, likelihood_values)
names(plot_data) <- c("x","y","z")
contour(plot_data, drawlabels = F)

optimal = optim(c(0.5,0.5), neg_likelihood)

dfll = data.frame(x = p_values, y = tau_values, llik = as.vector(likelihood_values))
dfll_dot = data.frame(x = optimal$par[1], y = optimal$par[2])

ggplot() + geom_contour(data = dfll, aes(x = x, y = y, z = llik, col = "Loglikelihood"), show.legend = T) +
  geom_point(data = dfll_dot, aes(x = x, y = y,col = "Computed optimal"), size = 4, alpha = 0.7) + labs(x = "p", y = "tau", title = "Loglikelihood of y, as function of p and tau")


# Exercise c)

p_est = optimal$par[1]
tau_est = optimal$par[1]
# First initiate normality constants and forward probabilities
norm_const = c(1/(dnorm(y[1],0,tau_est)*0.5+dnorm(y[1],1,tau_est)*0.5), rep(NA,249))
forward_prob = matrix(0,nrow=250,ncol = 2)
forward_prob[1,1] = norm_const[1]*0.5*dnorm(y[1],0,tau_est)
forward_prob[1,2] = norm_const[1]*0.5*dnorm(y[1],1,tau_est)
transition = matrix(NA, nrow = 250, ncol = 4)

# Iterate to find constants and probabilities
for (t in 2:250){
  norm_const[t] = 1/(dnorm(y[t],0,tau_est)*p_est*forward_prob[t-1,1] + (1-p_est)*dnorm(y[t],0,tau_est)*forward_prob[t-1,2] + dnorm(y[t],1,tau_est)*(1-p_est)*forward_prob[t-1,1] + dnorm(y[t],1,tau_est)*p_est*forward_prob[t-1,2])
  transition[t,1] = p_est*dnorm(y[t],0,tau_est)*forward_prob[t-1,1]*norm_const[t]
  transition[t,2] = (1-p_est)*dnorm(y[t],1,tau_est)*forward_prob[t-1,1]*norm_const[t]
  transition[t,3] = (1-p_est)*dnorm(y[t],0,tau_est)*forward_prob[t-1,2]*norm_const[t]
  transition[t,4] = p_est*dnorm(y[t],1,tau_est)*forward_prob[t-1,2]*norm_const[t]
  forward_prob[t,1] = transition[t,1] + transition[t,3]
  forward_prob[t,2] = transition[t,2] + transition[t,4]
}

# Initiate
backward_prob = matrix(0,nrow = 250, ncol = 2)
backward_prob[250,1] = forward_prob[250,1]
backward_prob[250,2] = forward_prob[250,2]
back_transition = matrix(NA, ncol = 4, nrow = 250)
propagation = matrix(NA, ncol = 4, nrow = 250)


for (t in 250:2){
  back_transition[t,1] = transition[t,1]/forward_prob[t,1]*backward_prob[t,1]
  back_transition[t,2] = transition[t,2]/forward_prob[t,2]*backward_prob[t,2]
  back_transition[t,3] = transition[t,3]/forward_prob[t,1]*backward_prob[t,1]
  back_transition[t,4] = transition[t,4]/forward_prob[t,2]*backward_prob[t,2]
  backward_prob[t-1,1] = back_transition[t,1] + back_transition[t,2]
  backward_prob[t-1,2] = back_transition[t,3] + back_transition[t,4]
  propagation[t,1] = back_transition[t,1]/backward_prob[t-1,1]
  propagation[t,2] = back_transition[t,2]/backward_prob[t-1,1]
  propagation[t,3] = back_transition[t,3]/backward_prob[t-1,2]
  propagation[t,4] = back_transition[t,4]/backward_prob[t-1,2]
}