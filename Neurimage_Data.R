## load the packages 
library(oro.nifti)
library(keras)
library(tidyverse)
library(plotly)
library(splines)
library(splines2)

## Read into image
img <- readANALYZE("/Users/chenhuan/Documents/桌面整理/签证/JHU/Research/EM_CurveRegression/EM_CurveRegression/tempFile.img")

index = which(img@.Data > 100, arr.ind = TRUE)
x = index[, 1] %>% as.vector
y = index[, 2] %>% as.vector
z = index[, 3] %>% as.vector

plot_ly() %>%
  add_trace(x = x, y = y, z = z, type = "scatter3d", mode = "markers", name = 'points', marker = list(size = 1, color = 'rgba(0, 0, 0, .9)', opacity = 0.4))

x = x / 128
y = y / 128
z = z / 128

###################################################################
################## The parameters to be tuned #####################
###################################################################
## Number of nodes 
N = length(x)
K = 200

## Degree of freedom for the splines
degree_free = 30

## Penalty coefficient
lambda1 = 1.5
lambda2 = 0

## The position of the curve
x_fix = c(0.35, 0.60)
y_fix = c(0.62, 0.57)
z_fix = c(0.14, 0.82)

##################################################################
##################################################################

## Basis for the spline
w <- seq(0, 2 * pi, by = 2 * pi / K)[1 : K]
B = cbind(1, bs(w, df = degree_free))
B_der = cbind(0, dbs(w, df = degree_free))
B_tilde = B[c(1,nrow(B)), ]


## PI is the matrix of p_ik
PI = matrix(1 / K, nrow = N, ncol = K)
PI_sum_old = rep(1 / K, K)

## sigma is the variance of the noise in the gaussian distribution
sigma_old = 1

## beta_x and beta_y are the coefficients for the splines to the x-axis and y-axis
beta_x_old = runif(degree_free + 1, -5, 5)
beta_y_old = runif(degree_free + 1, -5, 5)
beta_z_old = runif(degree_free + 1, -5, 5)

likelihood_store = c()

length_penalty = t(B_der) %*% B_der / K * 2 * pi 

## The procedure of the EM-Algorithm
for(t in 1 : 1000){
  ## items used during the EM procedure
  x.i.matrix = matrix(x,nrow=length(x),ncol=K,byrow=FALSE)
  x.k.matrix = matrix(B %*% beta_x_old, nrow = N, ncol = length(B %*% beta_x_old), byrow = TRUE)
  y.i.matrix = matrix(y,nrow = length(y), ncol = K, byrow = FALSE)
  y.k.matrix = matrix(B %*% beta_y_old, nrow = N, ncol = length(B %*% beta_y_old), byrow = TRUE)
  z.i.matrix = matrix(z, nrow = length(z), ncol = K, byrow = FALSE)
  z.k.matrix = matrix(B %*% beta_z_old, nrow = N, ncol = length(B %*% beta_z_old), byrow = TRUE)
  
  
  ## E-step
  PI = exp(-1 / as.numeric((2 * sigma_old)) * ((x.i.matrix - x.k.matrix) ^ 2 + (y.i.matrix - y.k.matrix)^ 2 + (z.i.matrix - z.k.matrix)^ 2)) %*% diag(PI_sum_old)
  PI = PI / apply(PI, 1, sum)

  
  ## M-step
  ## Update PI_sum
  PI_sum_new = 1 / N * apply(PI, 2, sum)
  
  ## Update sigma
  sigma_temp = 0
  sigma_temp = sum(((x.i.matrix - x.k.matrix)^2 + (y.i.matrix - y.k.matrix)^2 + (z.i.matrix - z.k.matrix)^2) * PI)
  sigma_new = 2 * sigma_temp / (3 * N)
  
  ## Update beta_x and beta_y
  B_XX = 0
  B_YY = 0
  B_ZZ = 0
  B_XY = 0
  
  for(i in 1 : N){
    #B_XY = B_XY + t(B) %*% diag(PI[i, ]) %*% B
    B_XX = B_XX + t(B) %*% as.matrix(PI[i, ]) * x[i] 
    B_YY = B_YY + t(B) %*% as.matrix(PI[i, ]) * y[i]
    B_ZZ = B_ZZ + t(B) %*% as.matrix(PI[i, ]) * z[i] 
  }
  
  diag_B = apply(PI, 2, sum) %>% diag
  B_XY = t(B) %*% diag_B %*% B
  #B_XX = apply(t(B) %*% (t(PI) %*% diag(x)), 1, sum)
  
  
  
  ## Inverse matrix for the estimation
  Inverse_M = solve(B_XY + lambda1 * length_penalty + lambda2 * diag(degree_free + 1))
  
  beta_x_new  = Inverse_M %*% B_XX
  beta_y_new  = Inverse_M %*% B_YY
  beta_z_new  = Inverse_M %*% B_ZZ
  
  ## Psu-inverse matrix for the estimation of coefficient
  
  Inverse_P = Inverse_M %*% t(B_tilde) %*% 
    solve(B_tilde %*% Inverse_M %*% t(B_tilde))
  
  beta_x_new = beta_x_new - 
    Inverse_P %*%
    (B_tilde %*% beta_x_new - x_fix)
  
  beta_y_new = beta_y_new - 
    Inverse_P %*%
    (B_tilde %*% beta_y_new - y_fix)
  
  beta_z_new = beta_z_new - 
    Inverse_P %*%
    (B_tilde %*% beta_z_new - z_fix)
  
  ## Computation of the log likelihood
  # likelihood = 0
  # for(i in 1 : N){
  #   likelihood_temp = 0
  #   for(k in 1 : K){
  #     likelihood_temp = likelihood_temp + PI_sum_new[k] * 1 / sigma_new * exp(-1/(2 * sigma_new) * ((x[i] - B[k, ] %*% beta_x_new )^2 + (y[i] - B[k, ] %*% beta_y_new)^2+ (z[i] - B[k, ] %*% beta_z_new)^2))
  #   }
  #   likelihood = likelihood +  log(likelihood_temp)
  # }

  PI_sum_old = PI_sum_new
  sigma_old = sigma_new
  beta_x_old = beta_x_new
  beta_y_old = beta_y_new
  beta_z_old = beta_z_new
  
   print(t)
  # likelihood_store = c(likelihood_store, likelihood)
  
  
  ## Save the runing result 
  #curve = list(X.fit = B %*% beta_x_new, Y.fit = B %*% beta_y_new, x = x, y = y)
  
  #saveRDS(curve, "res.rds")
  
  if(t %% 10 == 0){
    x.fit = B %*% beta_x_new
    y.fit = B %*% beta_y_new
    z.fit = B %*% beta_z_new
    plot_ly() %>%
      add_trace(x = x, y = y, z = z, type = "scatter3d", mode = "markers", name = 'points', marker = list(size = 1, color = 'rgba(0, 0, 0, .9)', opacity = 0.4)) %>%
      add_trace(x = as.vector(x.fit), y = as.vector(y.fit), z = as.vector(z.fit), type = "scatter3d", mode = "lines", name = "theoretical line", line = list(width = 2, color = 'rgba(255, 0, 0, .9)'))
    # par(mfrow = c(1,2))
    # scatter3D(x, y, z)
    # scatter3D(B %*% beta_x_new, B %*% beta_y_new, B %*% beta_z_new)
  }
  
}
dat = list(x = x, y = y)



