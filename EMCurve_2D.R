library(splines)
library(splines2)
library(tidyverse)
library(CVXR)
library(plotly)

## The penality of the regression coefficient 
lambda = 1

## N is the number of total points on the image
N = 400

## K is the number of total points on the curve 
K = 40

## The parameterization of the curve 
w <- seq(0, 2 * pi, by = 2 * pi / K)[1 : K]


## Generate of the Curve
X = 0
Y = 0
Z = 0
x.para = runif(1, -1, 1)
y.para = runif(1, -1, 1)
z.para = runif(1, -1, 1)

x = NULL
y = NULL
z = NULL


## Curve Generation with Spiral
coeff = runif(1, -2, 2)
X = X + x.para * w * sin(coeff * w)
Y = Y + y.para * w * cos(coeff * w)
Z = Z + z.para * w


# ## Generate the coefficients for the splines
# x.coeff = runif(degree_free, -1, 1)
# y.coeff = runif(degree_free, -1, 1)
# z.coeff = runif(degree_free, -1, 1)
# X = B %*% x.coeff
# Y = B %*% y.coeff
# Z = B %*% z.coeff


for(k in 1 : K){
  x = c(x, X[k] + rnorm(N / K, 0, 0.1))
  y = c(y, Y[k] + rnorm(N / K, 0, 0.1))
  z = c(z, Z[k] + rnorm(N / K, 0, 0.1))
}

x = x - min(0, min(x))
y = y - min(0, min(y))

x = x / max(abs(x))
y = y/ max(abs(y))
plot(x,y)

###################################################################
################## The parameters to be tuned #####################
###################################################################
## Number of nodes 
N = length(x)
K = 50

## Degree of freedom for the splines
degree_free = 20

## Penalty coefficient
lambda1 = 1
lambda2 = 0

## The position of the curve
x_fix = c(0.4, 0.01)
y_fix = c(0.4, 0.9)

##################################################################
##################################################################

## Apply EM Algorithm to the image
## Basis for the spline
w <- seq(0, 2 * pi, by = 2 * pi / K)[1 : K]
B = cbind(1, bs(w, df = degree_free))
B_der = cbind(0, dbs(w, df = degree_free))
B_tilde = B[c(1, nrow(B)), ]


## PI is the matrix of p_ik
PI = matrix(1 / K, nrow = N, ncol = K)
PI_sum_old = rep(1 / K, K)

## sigma is the variance of the noise in the gaussian distribution
sigma_old = 1

## beta_x and beta_y are the coefficients for the splines to the x-axis and y-axis
beta_x_old = runif(degree_free + 1, -5, 5)
beta_y_old = runif(degree_free + 1, -5, 5)

likelihood_store = c()

length_penalty = t(B_der) %*% B_der / K * 2 * pi 


## The procedure of the EM-Algorithm
for(t in 1 : 1000){
  ## E-step
  # for(i in 1 : N){
  #   for(k in 1 : K){
  #     PI[i,k] = exp(-1 / (2 * sigma_old) * ((x[i] - B[k, ] %*% beta_x_old) ^ 2 + (y[i] - B[k, ] %*% beta_y_old)^ 2)) * PI_sum_old[k]
  #   }
  #   PI = PI / apply(PI, 1, sum)
  # }
  
  
  
  ## items used during the EM procedure
  x.i.matrix = matrix(x,nrow=length(x),ncol=K,byrow=FALSE)
  x.k.matrix = matrix(B %*% beta_x_old, nrow = N, ncol = length(B %*% beta_x_old), byrow = TRUE)
  y.i.matrix = matrix(y,nrow = length(y), ncol = K, byrow = FALSE)
  y.k.matrix = matrix(B %*% beta_y_old, nrow = N, ncol = length(B %*% beta_y_old), byrow = TRUE)
  
  ## E-step
  PI = exp(-1 / as.numeric((2 * sigma_old)) * ((x.i.matrix - x.k.matrix) ^ 2 + (y.i.matrix - y.k.matrix)^ 2)) %*% diag(PI_sum_old)
  PI = PI / apply(PI, 1, sum)
  
  
  
  ## M-step
  ## Update PI_sum
  PI_sum_new = 1 / N * apply(PI, 2, sum)
  
  ## Update sigma
  sigma_temp = 0
  sigma_temp = sum(((x.i.matrix - x.k.matrix)^2 + (y.i.matrix - y.k.matrix)^2 ) * PI)
  sigma_new =  sigma_temp / (2 * N)
  
  
  ## Update beta_x and beta_y
  B_XX = 0
  B_YY = 0
  B_ZZ = 0
  B_XY = 0
  
  for(i in 1 : N){
    #B_XY = B_XY + t(B) %*% diag(PI[i, ]) %*% B
    B_XX = B_XX + t(B) %*% as.matrix(PI[i, ]) * x[i] 
    B_YY = B_YY + t(B) %*% as.matrix(PI[i, ]) * y[i]
  }
  
  diag_B = apply(PI, 2, sum) %>% diag
  B_XY = t(B) %*% diag_B %*% B
  #B_XX = apply(t(B) %*% (t(PI) %*% diag(x)), 1, sum)
  
  
  
  ## Inverse matrix for the estimation
  Inverse_M = solve(B_XY + lambda1 * length_penalty + lambda2 * diag(degree_free + 1))
  
  beta_x_new  = Inverse_M %*% B_XX
  beta_y_new  = Inverse_M %*% B_YY
  
  ## Psu-inverse matrix for the estimation of coefficient
  
  Inverse_P = Inverse_M %*% t(B_tilde) %*% 
    solve(B_tilde %*% Inverse_M %*% t(B_tilde))
  
  beta_x_new = beta_x_new - 
    Inverse_P %*%
    (B_tilde %*% beta_x_new - x_fix)
  
  beta_y_new = beta_y_new - 
    Inverse_P %*%
    (B_tilde %*% beta_y_new - y_fix)
  
  # ## CVXR method to do the estimation
  # beta_x_cvx = Variable(degree_free + 1)
  # beta_y_cvx = Variable(degree_free + 1)
  # 
  # loss_x1 = sum(x * x)
  # loss_y1 = sum(y * y)
  # loss_x2 = 0
  # loss_y2 = 0
  # loss_x3 = 0
  # loss_y3 = 0
  # for(i in 1 : N){
  #   for(k in 1 : K){
  #     loss_x2 = loss_x2 + B[k, ] * PI[i, k] * x[i]
  #     loss_y2 = loss_y2 + B[k, ] * PI[i, k] * y[i]
  #     loss_x3 = loss_x3 + B[k, ] %*% t(B[k, ]) * PI[i, k]
  #     loss_y3 = loss_x3
  #   }
  # }
  # loss_penalty_x = 0
  # loss_penalty_y = 0
  # for(k in 1 : K){
  #   loss_penalty_x = loss_penalty_x + B_der[k, ] %*% t(B_der[k, ]) 
  #   loss_penalty_y = loss_penalty_y + B_der[k, ] %*% t(B_der[k, ]) 
  # }
  # loss_penalty = (quad_form(beta_x_cvx, loss_penalty_x) + quad_form(beta_y_cvx, loss_penalty_y)) * 2 * pi / K
  # 
  # loss =
  #   loss_x1  - 2 *  t(beta_x_cvx) %*% loss_x2 + quad_form(beta_x_cvx, loss_x3) +
  #   loss_y1  - 2 *  t(beta_y_cvx) %*% loss_y2 + quad_form(beta_y_cvx, loss_y3) + 
  #   loss_penalty
  # 
  # objective <- Minimize(loss)
  # prob <- Problem(objective)
  # CVXR_result <- solve(prob)
  # 
  # beta_x_new = CVXR_result$getValue(beta_x_cvx)
  # beta_y_new = CVXR_result$getValue(beta_y_cvx)
  
  
  ## Computation of the log likelihood
  likelihood = 0
  for(i in 1 : N){
    likelihood_temp = 0
    for(k in 1 : K){
      likelihood_temp = likelihood_temp + PI_sum_new[k] * 1 / sigma_new * exp(-1/(2 * sigma_new) * ((x[i] - B[k, ] %*% beta_x_new )^2 + (y[i] - B[k, ] %*% beta_y_new)^2))
    }
    likelihood = likelihood +  log(likelihood_temp)
  }
  
  PI_sum_old = PI_sum_new
  sigma_old = sigma_new
  beta_x_old = beta_x_new
  beta_y_old = beta_y_new
  
  print(likelihood)
  likelihood_store = c(likelihood_store, likelihood)
  
  
  ## Save the runing result 
  #curve = list(X.fit = B %*% beta_x_new, Y.fit = B %*% beta_y_new, x = x, y = y)
  
  if(t %% 10 == 0){
    plot(x,y, xlim = c(0,1), ylim = c(0,1))
    lines(B %*% beta_x_new, B %*% beta_y_new, type = "l", col = "red")
  }
  
}
dat = list(x = x, y = y)
#saveRDS(dat, "data.rds")



plot(x,y, xlim = c(0,1), ylim = c(0,1))
lines(B %*% beta_x_new, B %*% beta_y_new, type = "l", col = "red")





