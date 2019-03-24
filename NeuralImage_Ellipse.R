## load the packages 
library(oro.nifti)
library(keras)
library(tidyverse)
library(plotly)
library(splines)
library(splines2)
library(pracma)


## Read into image
img <- readANALYZE("/Users/chenhuan/Documents/桌面整理/签证/JHU/Research/EM_CurveRegression/EM_CurveRegression/tempFile.img")

index = which(img@.Data > 150, arr.ind = TRUE)
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
K = 50

## Degree of freedom for the splines
degree_free = 20

## Penalty coefficient
lambda1 = 40
lambda2 = 0

## The position of the curve
#x_fix = c(0.35, 0.60)
#y_fix = c(0.62, 0.57)
#z_fix = c(0.14, 0.82)

# x_fix = c(0.35, 0.57)
# y_fix = c(0.60, 0.54)
# z_fix = c(0.17, 0.83)

x_fix = c(0.35, 0.57)
y_fix = c(0.60, 0.57)
z_fix = c(0.17, 0.85)


alpha = 2
beta = 0.0001

T = 50

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
sigma_old_x = 1
sigma_old_y = 1
sigma_old_z = 1

## beta_x and beta_y are the coefficients for the splines to the x-axis and y-axis
beta_x_old = runif(degree_free + 1, -5, 5) * 0
beta_y_old = runif(degree_free + 1, -5, 5) * 0
beta_z_old = runif(degree_free + 1, -5, 5) * 0

likelihood_store = c()

length_penalty = t(B_der) %*% B_der / K * 2 * pi 

## The procedure of the EM-Algorithm
for(t in 1 : T){
  
  ## items used during the EM procedure
  x.i.matrix = matrix(x,nrow=length(x),ncol=K,byrow=FALSE)
  x.k.matrix = matrix(B %*% beta_x_old, nrow = N, ncol = length(B %*% beta_x_old), byrow = TRUE)
  y.i.matrix = matrix(y,nrow = length(y), ncol = K, byrow = FALSE)
  y.k.matrix = matrix(B %*% beta_y_old, nrow = N, ncol = length(B %*% beta_y_old), byrow = TRUE)
  z.i.matrix = matrix(z, nrow = length(z), ncol = K, byrow = FALSE)
  z.k.matrix = matrix(B %*% beta_z_old, nrow = N, ncol = length(B %*% beta_z_old), byrow = TRUE)
  
  
  ## E-step
  sigma_matrix_x = matrix(sigma_old_x, nrow = N, ncol = K, byrow = TRUE)
  sigma_matrix_y = matrix(sigma_old_y, nrow = N, ncol = K, byrow = TRUE)
  sigma_matrix_z = matrix(sigma_old_z, nrow = N, ncol = K, byrow = TRUE)
  
  ## Compute the conditional expectation PI
  PI = (1 / sqrt(sigma_matrix_x * sigma_matrix_y * sigma_matrix_z) *  exp(-1 / (2 * sigma_matrix_x) * (x.i.matrix - x.k.matrix) ^ 2 - 1 / (2 * sigma_matrix_y) * (y.i.matrix - y.k.matrix)^ 2 - 1 / (2 * sigma_matrix_z) * (z.i.matrix - z.k.matrix)^ 2)) %*% diag(PI_sum_old)
  PI = PI / apply(PI, 1, sum)
  
  # PI_log =  1 -3 / 2 * log(sigma_matrix) - 1 / (2 * sigma_matrix) * ((x.i.matrix - x.k.matrix) ^ 2 + (y.i.matrix - y.k.matrix)^ 2 + (z.i.matrix - z.k.matrix)^ 2)  +  matrix(log(PI_sum_old), nrow = N, ncol = length(B %*% beta_z_old), byrow = TRUE)
  # 
  # PI = exp(PI_log) / apply(exp(PI_log), 1, sum)
  
  
  ## M-step
  ## Update PI_sum
  PI_sum_new = 1 / N * apply(PI, 2, sum)
  
  ## Update sigma
  sigma_temp_x = apply(((x.i.matrix - x.k.matrix)^2 ) * PI, 2, sum)
  sigma_new_x = (sigma_temp_x) / (apply(PI, 2, sum))
  
  sigma_temp_y = apply(((y.i.matrix - y.k.matrix)^2 ) * PI, 2, sum)
  sigma_new_y = (sigma_temp_y) / (apply(PI, 2, sum))
  
  sigma_temp_z = apply(((z.i.matrix - z.k.matrix)^2 ) * PI, 2, sum)
  sigma_new_z = (sigma_temp_z) / (apply(PI, 2, sum))
  # sigma_new = exp(log(sigma_temp) - log(3) - log(apply(PI, 2, sum)))
  
  ## PI matrix after adjusting for the variance component
  #PI_reg = PI / matrix(sigma_new, nrow = N, ncol = K, byrow = TRUE)
  PI_reg = PI
  ## Update beta_x and beta_y
  B_XX = 0
  B_YY = 0
  B_ZZ = 0
  B_XY = 0
  
  for(i in 1 : N){
    #B_XY = B_XY + t(B) %*% diag(PI[i, ]) %*% B
    B_XX = B_XX + t(B) %*% as.matrix(PI_reg[i, ]) * x[i] 
    B_YY = B_YY + t(B) %*% as.matrix(PI_reg[i, ]) * y[i]
    B_ZZ = B_ZZ + t(B) %*% as.matrix(PI_reg[i, ]) * z[i] 
  }
  
  diag_B = apply(PI_reg, 2, sum) %>% diag
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
  sigma_old_x = sigma_new_x
  sigma_old_y = sigma_new_y
  sigma_old_z = sigma_new_z
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
      add_trace(x = as.vector(x.fit), y = as.vector(y.fit), z = as.vector(z.fit), type = "scatter3d", mode = "lines", name = "theoretical line", line = list(width = 5, color = 'rgba(255, 0, 0, .9)'))
    # par(mfrow = c(1,2))
    # scatter3D(x, y, z)
    # scatter3D(B %*% beta_x_new, B %*% beta_y_new, B %*% beta_z_new)
    #print(x.fit)
    
  }
  print(sigma_new_x)
  
}

## Plot the circle result
plt_store = 
  plot_ly() %>%
  add_trace(x = x, y = y, z = z, type = "scatter3d", mode = "markers", name = 'points', marker = list(size = 1, color = 'rgba(0, 0, 0, .9)', opacity = 0.4)) %>%
  add_trace(x = as.vector(x.fit), y = as.vector(y.fit), z = as.vector(z.fit), type = "scatter3d", mode = "lines", name = "theoretical line", line = list(width = 5, color = 'rgba(255, 0, 0, .9)'))

for(i in 1 : 500){
  m = Circle_Plot(c(x.fit[i],y.fit[i],z.fit[i]), diag(c(sigma_new_x[i], sigma_new_y[i], sigma_new_z[i])), 2, c(B_der[i, ] %*%beta_x_new,B_der[i, ] %*%beta_y_new, B_der[i, ] %*%beta_z_new))
  plt_store = plt_store %>% add_trace(x = m[1, ], y = m[2, ], z = m[3, ], type = "scatter3d", mode = "lines", name = "theoretical line", line = list(width = 5, color = paste0('rgba(',  0.5 * i,', 0, 0, .9)')))
  
}
plt_store %>% layout(showlegend = FALSE)


col_value = sigma_new_x + sigma_new_y + sigma_new_z
index_value = sort(col_value)
match(index_value, col_value)


