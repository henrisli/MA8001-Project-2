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

ggplot(df, aes(x = 1:250, y = x)) + geom_line() + labs(x = "i", y = "x", title = "Generated sample for x") + theme_classic(base_size = 19)

ggplot(df, aes(x = 1:250, y = y)) + geom_point() + labs(x = "i", y = "y", title = "Generated sample for y") + theme_classic(base_size = 19)

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
  geom_point(data = dfll_dot, aes(x = x, y = y,col = "Computed optimal"), size = 4, alpha = 0.7) + labs(x = "p", y = "tau", title = "Loglikelihood of y, as function of p and tau") + theme_classic()


# Exercise c)

p_est = optimal$par[1]
tau_est = optimal$par[2]
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

df2 = data.frame(prob_0 = backward_prob[,1], prob_1 = backward_prob[,2])
ggplot(df2, aes(x = 1:250, y = prob_1)) + geom_line() + labs(x = "i", y = "prob x_i = 1", title = "Computed marginal probabilities for x_i = 1") + theme_classic(base_size = 19)

x_est = rep(NA, 250)
x_est[1] = ifelse(runif(1)<backward_prob[1,2], 1, 0)
for (t in 2:250){
  x_est[t] = ifelse(runif(1) < propagation[t,as.integer(2*(1+x_est[t-1]))], 1, 0)
}

df3 = data.frame(x = x_est)

ggplot(df3, aes(x = 1:250, y = x)) + geom_line() + labs(x = "i", y = "x", title = "Simulated values for x based on FB-algorithm") + theme_classic(base_size = 19)


# Exercise d)

df4 = data.frame(markov = round(backward_prob[,2]), indep = ifelse(y<0.5, 0, 1))
ggplot(df4, aes(x = 1:250)) + geom_line(aes(y = markov)) + labs(x = "i", y = "x", title = "Predicted values for x based on Markov property") + theme_classic(base_size = 19)
ggplot(df4, aes(x = 1:250)) + geom_line(aes(y = indep)) + labs(x = "x", title = "Predicted values for x based on independence ") + theme_classic(base_size = 19)


df2 = data.frame(prob_0 = backward_prob[,1], prob_1 = backward_prob[,2], prob_1_indep = dnorm(y,1,tau)/(dnorm(y,1,tau)+dnorm(y,0,tau)))
ggplot(df2, aes(x = 1:250)) + geom_line(size = 1.2, aes(y = prob_1, col = "Markov")) + geom_line(size = 0.1, aes(y = prob_1_indep, col = "Independent")) + labs(x = "i", y = "prob x_i = 1", title = "Probability of x_i = 1 for Markov and independence assumptions") + theme_classic()
