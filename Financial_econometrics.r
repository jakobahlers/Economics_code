#####################################################################
##########   Financial Econometrics Exam - 2021/2020    #############
##########   Flow number: 55                            #############
#####################################################################
#####################################################################
library(rugarch)
library(moments) # loading package to compute moments
library("rmgarch") # loading rmgarch package 


#############################################
#########   Computational part  #############
#############################################
#############################################
##########      Question 1      #############
#############################################

tse_tsu <- function(i_T, i_P, i_M, d_omega, d_alpha, d_beta, d_a, d_b, m_R) {

    # This function's main purpose is to simulate observations from the Tse 
    # and Tsui (2002) model. 
    # It takes the following inputs: 
        # i_T; periods
        # i_P; number of assets. Be aware that this function is coded, 
             # so it only can take P=2. This is off course a limitation.
             # This could be generalized, so the simulation could deal 
             # with P>2.
        # i_M; Up till period M, we have that R_t = R
        # Parameters of the model; d_omega, d_alpha, d_beta, d_a, d_b
        # m_R; Correlation matrix up till and including period M. 
    # The function outputs:
        # m_Y; matrix with simulated returns
        # a_Sigma; array with covariance matrix for each time period
        # a_R; time-varying correlation matrix
        # a_Psi; "driving force" - see slide 37, lecture 10, for further description
        # a_D; diagonal matrix with generic element sigma_i,t --> See slide 30, lecture 10
        # m_sigma2; Matrix with volatility process for each asset
    # Using "expm" package which enables us to take the squareroot of Sigma, when computing y. 
    require("expm")

    # Storage
    m_Y <- matrix(data = NA, nrow = i_T, ncol = i_P, byrow = TRUE)
    a_D <- array(0, dim = c(i_P, i_P, i_T))
    a_R <- array(0, dim = c(i_P, i_P, i_T))
    a_Sigma <- array(0, dim = c(i_P, i_P, i_T))
    m_sigma2 <- matrix(data = NA, nrow = i_T, ncol = i_P, byrow = TRUE)
    a_Psi <- array(0, dim = c(i_P, i_P, i_T))

    # Creating epsilon
    m_eps <- matrix(data = rnorm(i_T * i_P), nrow = i_T, ncol = i_P, byrow = TRUE)

    # Initializing D by unconditional value of sigma_1^2 = omega / (1- alpha - beta)
    m_sigma2[1, ] <- d_omega / (1.0 - d_alpha - d_beta)
    # Inserting elements into the diagonal of D matrix. 
    diag(a_D[,,1]) <- sqrt(m_sigma2[1, ])

    # Initializing R
    a_R[,,1] <- m_R

    # Computing Sigma_1
    a_Sigma[,,1] <- a_D[,,1] %*% a_R[,,1] %*% a_D[,,1]

    # Initializing Y
    m_Y[1, ] <- sqrtm(a_Sigma[,,1]) %*% m_eps[1,]

    # Loop for time up till period M, where R_t = R
    for (i_t in seq(2,i_M)) {
        # Computing sigma2 and inserting into sqrt(sigma2) in the diagonal of D.
        m_sigma2[i_t, ] <- d_omega + d_alpha * m_Y[i_t - 1, ]^2 + d_beta * m_sigma2[i_t - 1, ]
        diag(a_D[,,i_t]) <- sqrt(m_sigma2[i_t, ])

        # Correlation matrix, which is constant for t<=M
        a_R[,,i_t] <- m_R

        # Computing Sigma
        a_Sigma[,,i_t] <- a_D[,,i_t] %*% a_R[,,i_t] %*% a_D[,,i_t]

        # Computing Y
        m_Y[i_t, ] <- sqrtm(a_Sigma[,,i_t]) %*% m_eps[i_t,]

    }

    # Loop for time M up till period T, where R_t is time varying
    for (i_t in seq(i_M + 1, i_T)) {
        # Computing sigma2 and inserting into sqrt(sigma2) in the diagonal of D.
        m_sigma2[i_t, ] <- d_omega + d_alpha * m_Y[i_t - 1, ]^2 + d_beta * m_sigma2[i_t - 1, ]
        diag(a_D[,,i_t]) <- sqrt(m_sigma2[i_t, ])

        # Computing eta
        m_eta <- m_Y / sqrt(m_sigma2)

        # Computing Psi
        # For the off-diagonal where "i" is not equal to "j"
        nominator <- sum(m_eta[seq(i_t - 1 ,i_t - i_M) , 1] * m_eta[seq(i_t - 1 ,i_t - i_M), 2])
        denominator <- sqrt(sum(m_eta[seq(i_t - 1 ,i_t - i_M) , 1]^2) * sum(m_eta[seq(i_t - 1 ,i_t - i_M) , 2]^2))
        d_psi <- nominator / denominator
        # Setting into Psi
        a_Psi[1, 2, i_t - 1] <- d_psi
        a_Psi[2, 1, i_t - 1] <- d_psi

        # Note that it is not nessesary to compute this for i = j, 
        # which can be seen from the formula for phi. 
        # Therefore the code would be more efficient if we just put 1.0 
        # into the diagonal of Psi instead of computing the following. 
        # For the diagonal where i = j (asset 1)
        nominator <- sum(m_eta[seq(i_t - 1 ,i_t - i_M) , 1] * m_eta[seq(i_t - 1 ,i_t - i_M), 1])
        denominator <- sqrt(sum(m_eta[seq(i_t - 1 ,i_t - i_M) , 1]^2) * sum(m_eta[seq(i_t - 1 ,i_t - i_M) , 1]^2))
        d_psi <- nominator / denominator
        a_Psi[1, 1, i_t - 1] <- d_psi

        # and again i = j (asset 2)
        nominator <- sum(m_eta[seq(i_t - 1 ,i_t - i_M) , 2] * m_eta[seq(i_t - 1 ,i_t - i_M), 2])
        denominator <- sqrt(sum(m_eta[seq(i_t - 1 ,i_t - i_M) , 2]^2) * sum(m_eta[seq(i_t - 1 ,i_t - i_M) , 2]^2))
        d_psi <- nominator / denominator
        a_Psi[2, 2, i_t - 1] <- d_psi

        # Updating R according to update equation
        a_R[,,i_t] <- (1 - d_a - d_b) * m_R + d_a * a_Psi[,, i_t - 1] + d_b * a_R[,,i_t - 1]

        # Computing Sigma
        a_Sigma[,,i_t] <- a_D[,,i_t] %*% a_R[,,i_t] %*% a_D[,,i_t]

        # Computing Y
        m_Y[i_t, ] <- sqrtm(a_Sigma[,,i_t]) %*% m_eps[i_t,]

    }
    
    # Creating list which will be output of function
    l_out <- list()
    l_out[["m_Y"]] <- m_Y
    l_out[["a_Sigma"]] <- a_Sigma
    l_out[["a_R"]] <- a_R
    l_out[["a_Psi"]] <- a_Psi
    l_out[["a_D"]] <- a_D
    l_out[["m_sigma2"]] <- m_sigma2
    l_out[["m_eta"]] <- m_eta
    return(l_out)

}



#############################################
##########      Question 2      #############
#############################################

# Note that the following functions follows a similar procedure as 
# from exercise set 6. 
tse_tsui_Filter <- function(m_Eta, d_a, d_b, m_Q, d_M = 5) {

    # Main purpose of this is to compute correlation component of log-likelihood 
    # and correlation matrices for each time period t. 
    # The function is used in estimation function in order to estimate Tse and
    # Tsui (2002) model. 
    # Inputs: 
        # m_Eta: Matrix with elements: y_j,t / sigma_j,t for j=1,...,p
        # d_a and d_b: Present in updating equation for R (correlation matrix)
        # m_Q: used for initialization of R
        # d_M: number of periods, where R is constant.

    i_P = ncol(m_Eta)
    i_T = nrow(m_Eta)

    # initialize storage arrays and matrices
    a_R = array(0, dim = c(i_P, i_P, i_T))
    a_D <- array(0, dim = c(i_P, i_P, i_T))
    a_Psi <- array(0, dim = c(i_P, i_P, i_T))

    ## initialization at the unconditional cor
    m_R <- m_Q
    a_R[,, 1] = m_R

    # Initialize LLK
    d_LLK <- t(m_Eta[1, ]) %*% solve(a_R[,, 1]) %*% m_Eta[1, ] - t(m_Eta[1, ]) %*% m_Eta[1, ] + log(det(a_R[,, 1]))
    

    for (i_t in 2:d_M) {

        # R matrix, which is constant up to and including period M
        a_R[,,i_t] <- m_R

        # augment the likelihood value (correlation component, slide 35, lecture 10)
        d_LLK <- d_LLK + t(m_Eta[i_t, ]) %*% solve(a_R[,, i_t]) %*% m_Eta[i_t, ] - t(m_Eta[i_t, ]) %*% m_Eta[i_t, ] + log(det(a_R[,, i_t]))
    }


   # Loop for time M up till period T, where R_t is time varying
    for (i_t in seq(i_M + 1, i_T)) {

        # Computing Psi
        # Note that this code could be optimized if the following were written vectorized. 
        # Further it should be noticed that the function do not work if P>2. 
        # In the empirical part of the code this function is generalized, 
        # so it works for P>2. 
        # For the off-diagonal where "i" is not equal to "j"
        nominator <- sum(m_Eta[seq(i_t - 1 ,i_t - i_M) , 1] * m_Eta[seq(i_t - 1 ,i_t - i_M), 2])
        denominator <- sqrt(sum(m_Eta[seq(i_t - 1 ,i_t - i_M) , 1]^2) * sum(m_Eta[seq(i_t - 1 ,i_t - i_M) , 2]^2))
        d_psi <- nominator / denominator
        # Setting into Psi
        a_Psi[1, 2, i_t - 1] <- d_psi
        a_Psi[2, 1, i_t - 1] <- d_psi

        # Note that here I use the knowledge that we have 1 on the diagonals of Psi
        # For the diagonal where i = j (asset 1)
        a_Psi[1, 1, i_t - 1] <- 1

        # and again i = j (asset 2)
        a_Psi[2, 2, i_t - 1] <- 1

        # Updating R
        a_R[,,i_t] <- (1 - d_a - d_b) * m_R + d_a * a_Psi[,, i_t - 1] + d_b * a_R[,,i_t - 1]

        # augment the likelihood value (correlation component, slide 35, lecture 10)
        d_LLK <- d_LLK + t(m_Eta[i_t, ]) %*% solve(a_R[,, i_t]) %*% m_Eta[i_t, ] - t(m_Eta[i_t, ]) %*% m_Eta[i_t, ] + log(det(a_R[,, i_t]))

    }

    # Creating list which should be returned by function
    l_Out = list()
    # Including the -1/2 term in the likelihood evaluation
    #see the equations in the corresponding lecture
    l_Out[["d_LLK"]] = -0.5 * d_LLK
    l_Out[["a_R"]] = a_R

    return(l_Out)

}


# Objective function which we would like to maximize. 
objective <- function(m_Eta, v_Par, m_Q) {
   
    Filter = tse_tsui_Filter(m_Eta, v_Par[1], v_Par[2], m_Q)
    # Computing negative log-likelihood
    d_NLLK = -as.numeric(Filter$d_LLK)
    return(d_NLLK)

}


Estimate_tse_tsui <- function(m_Y, d_M = 5) {

    # The main purpose of this function is to estimate model of Tse and Tsui (2002). 
    # Input: 
        # m_Y: Series of return
        # d_M: We have that R is constant up to and including time-period M
    # It outputs a list with several results. 
        # d_LLK: Log-likelihood
        # m_Coef: Estimated alpha and beta
        # v_Par: Estimated a and b, which is present in the updating equation for R. 
        # m_Sigma: Volatilities
        # a_R: Correlation matrices for each t
        # m_Eta: Residuals
        # m_Z
        # BIC: Bayesian information criteria

    ## estimate the marginal GARCH models
    require(rugarch) #loading the rugarch package
    require(Rsolnp)  #loading the Rsolnp package

    #Marginal garch specifications
    SpecGARCH = ugarchspec(mean.model = list(armaOrder = c(0, 0)))
    m_spec <- multispec( replicate(ncol(m_Y), SpecGARCH) )

    # Fitting garch models
    l_Fit_univariate <- multifit(multispec = m_spec, data = m_Y)

    # Computing residuals
    m_Eta <- residuals(l_Fit_univariate, standardize = TRUE)

    # Initial parameters
    v_Par = c(0.003, 0.995)

    # Setting empirical correlation from data
    m_Q <- cor(m_Y)

    # Maximizing likelihood
    optimizer <- solnp(v_Par, objective, 
                        ineqfun = function(v_Par, ...) {
                        sum(v_Par)
                        },
                        ineqLB = 1e-4, ineqUB = 0.999,
                        LB = c(1e-4, 1e-4), UB = c(0.999, 0.999),
                        m_Eta = m_Eta, m_Q = m_Q)

    # Extract parameters
    v_Par = optimizer$pars

    #Filter the dynamic correlation using the estimated parameters
    Filter = tse_tsui_Filter(m_Eta, v_Par[1], v_Par[2], m_Q, d_M)

    #extract univariate volatilities
    m_Sigma = sigma(l_Fit_univariate)

    #extract univariate estimated parameters
    m_Coef = coef(l_Fit_univariate)

    #compute the likelihood of the volatility  part
    d_LLK_V = likelihood(l_Fit_univariate)

    #compute the total likelihood
    a_R = Filter[["a_R"]]
    a_Cov = array(NA, dim = dim(a_R))
    i_T = dim(a_R)[3]
    i_P = dim(a_R)[2]
    vol_LLK = 0

    for (t in 1:i_T) {

    a_Cov[,,t] = diag(m_Sigma[t, ]) %*% a_R[,,t] %*% diag(m_Sigma[t, ])
    vol_LLK <- vol_LLK + i_P * log(2 * pi) + log(det(a_Cov[,, t])) + as.numeric(m_Y[t, ]) %*% solve(a_Cov[,,t]) %*% as.numeric(t(m_Y[t, ]))
   
    }

    # Log-likelihood correlation part
    LLK_cor <- Filter[["d_LLK"]]

    # Total log likelihood likelihood
    d_LLK <- -0.5 * (LLK_cor + vol_LLK)

    ## Compute z_t
    i_T = nrow(m_Y)
    # Storage
    m_Z = matrix(0, i_T, ncol(m_Y))

    for (t in 1:i_T) {
    m_Z[t, ] = diag(1/m_Sigma[t, ]) %*% solve(chol(a_R[,,t])) %*% as.numeric(m_Y[t, ])
    }

    # Computing BIC (Bayesian information criteria)
    BIC = log(i_T) * 8 - 2 * d_LLK

    # Creating list which should be returned
    l_Out = list()

    #output the results
    l_Out[["d_LLK"]] = d_LLK
    l_Out[["m_Coef"]] = m_Coef
    l_Out[["v_Par"]] = v_Par
    l_Out[["m_Sigma"]] = m_Sigma
    l_Out[["a_R"]] = a_R
    l_Out[["m_Eta"]] = m_Eta
    l_Out[["m_Z"]] = m_Z
    l_Out[["BIC"]] = BIC

    return(l_Out)

}

# As pointed out in the code it should be noticed that this code above do NOT work properly, 
# if P>2. This is not a problem when solving the computational part of the exam. 
# It though becomes a problem in the empirical part, where the exercise asks to fit the 
# Tse and Tsui (2002) model on 3 series of returns from the Dow Jones 30 index. In the empirical 
# part I have rewritten the code such that it works for P>2. This is done by rewriting the
# computation of Psi. 
# The reason for why I have chosen not to rewrite the code above for this question 2 
# is that it runs faster compared to the estimation function written in the empirical part. 
# This is probably due to the fact that when vectorizing the computation of Psi, the function
# needs to invert several matrices in a loop which runs for all time periods of consideration. 


#############################################
##########      Question 3     #############
#############################################

############
#### a) ####
############
# Simulating with T=1000 from the function created in question 1. 
# Setting parameters
d_omega <- 0.01
d_alpha <- 0.04
d_beta <- 0.95
d_a <- 0.003
d_b <- 0.995
i_M <- 5
i_P <- 2
i_T <- 1000
m_R <- matrix(data = c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2, byrow = TRUE)

# Setting the seed to 123
set.seed(123)
sim_data <- tse_tsu(i_T, i_P, i_M, d_omega, d_alpha, d_beta, d_a, d_b, m_R)

# I define variables we wish to plot
m_Y_1 <- sim_data$m_Y[ , 1]
m_Y_2 <- sim_data$m_Y[ , 2]
v_var_1 <- sim_data$a_Sigma[1,1,]
v_var_2 <- sim_data$a_Sigma[2,2,]
v_cor <- sim_data$a_R[1,2,]

# Plotting series of simulated returns, time-varying variances
# and time-varying correlations. 
par(mfrow=c(3,2))
plot.ts(m_Y_1, main = "Returns - Asset 1")
plot.ts(m_Y_2, main = "Returns - Asset 2")
plot.ts(v_var_1, main = "Time varying variance - Asset 1")
plot.ts(v_var_2, main = "Time varying variance - Asset 2")
plot.ts(v_cor, main = "Time varying correlation")



############
#### b) ####
############
# Estimating parameters of the model using function from point 2). 
# First I define Y
m_Y <- sim_data$m_Y
# Estimating model
tse_tsu_est <- Estimate_tse_tsui(m_Y)
# Pulling out estimates
m_estimates <- tse_tsu_est$m_Coef
# Observing for asset 1: 
    # alpha = 0.03345539
    # beta = 0.92667594
# Observing for asset 2:
    # alpha = 0.030412685
    # beta = 0.960402395
v_par_est <- tse_tsu_est$v_Par
# Observing for a and b: 
    # a = 0.009911028
    # b = 0.974471692

# Remembering that the true parameters for a and b is set to: 
    # a = 0.003
    # b = 0.995
v_true_parameters <- c(0.003, 0.995)

# Setting results up in a matrix, which makes it easier to compare estimated and true parameters
m_param_compare <- matrix(data = c(v_par_est, v_true_parameters), nrow = 2, ncol = 2, byrow = FALSE)
colnames(m_param_compare) <- c("Estimated", "TRUE")
rownames(m_param_compare) <- c("a", "b") 
#     Estimated  TRUE
# a 0.009911028 0.003
# b 0.974471692 0.995



#############################################
#########     Empirical part    #############
#############################################
#############################################
##########      Question 1      #############
#############################################
# Loading data from dji30ret
data("dji30ret")
colnames(dji30ret) # printing tickers
# Defining matrix with data
# I choose the first 3 series from dji30ret and multiply by 100
m_AA <- dji30ret[1:500, c("AA")] * 100
m_AXP <- dji30ret[1:500, c("AXP")] * 100
m_BA <- dji30ret[1:500, c("BA")] * 100

# Removing mean from series
m_AA <- m_AA - mean(m_AA)
m_AXP <- m_AXP - mean(m_AXP)
m_BA <- m_BA - mean(m_BA)

# Combining returns into matrix
m_returns <- matrix(data = c(m_AA, m_AXP, m_BA), nrow = length(m_AA), ncol = 3, byrow = FALSE)
colnames(m_returns) = c("AA", "AXP", "BA") # setting names for colunms

# Descripitve statistics of the three series
# First we compute correlation matrix for the tree series
m_correlation <- cor(m_returns) 
# Observing following correlations:
    # AA and AXP: 0.550447
    # AA and BA: 0.409803
    # AXP and BA: 0.5559013
# Observing that the correlation between AA/AXP and AXP/BA are almost similar
# around 0.55, while the correlation between AA/BA is a bit smaller.  

# Considering variance for each stock
m_variance <- var(m_returns)
# We get the variance for each stock on the diagonal of the matrix; m_variance
    # AA: 6.241963
    # AXP: 7.831274
    # BA: 3.200196
# When comparing the variance of the different returns series, it is 
# clear that AA and AXP have a greater variance compared to BA. 
# This greater variance of AA and AXP could imply that these return
# series are more volatile in comparison with the return series of BA. 

# Minimum and maximum of each stock
m_mm = c(max(m_AA), min(m_AA), max(m_AXP), min(m_AXP), max(m_BA), min(m_BA))
# Collecting into matrix
m_min_max <- matrix(data = m_mm, nrow = 2, ncol = 3, byrow = FALSE)
colnames(m_min_max) = c("AA", "AXP", "BA") # setting names for colunms
rownames(m_min_max) = c("max", "min") # setting names for rows
# We observe
#             AA       AXP        BA
# max   7.681382  17.15468  14.16899
# min -27.535632 -30.30871 -12.48065
# These observations fits nicely with the variance measure from previosly.
# We observe more extreme maximum and minimum return realizations for AA and AXP
# in comparison with the BA return series. 

# Next I will consider skewness and kurtosis of the different return series. 
m_skew_kurt = c(skewness(m_AA), kurtosis(m_AA), skewness(m_AXP), kurtosis(m_AXP), skewness(m_BA), kurtosis(m_BA))
m_momemts <- matrix(data = m_skew_kurt, nrow = 2, ncol = 3, byrow = FALSE)
colnames(m_momemts) = c("AA", "AXP", "BA") # setting names for colunms
rownames(m_momemts) = c("skewness", "kurtosis") # setting row names
#                AA       AXP         BA
# skewness -3.03780 -2.629072  0.2119338
# kurtosis 34.13447 35.091222 15.7883121
# We see that generally the skewness for the return series are negative, 
# which is typical for financial returns. It is though noticed that the skewness
# for BA is slightly positive. When extending the time horizon under consideration
# by another 500 periods this skewness though turns negative. 
# Kurtosis is generally large for all time series, but though largest for AA. 
# Note that kurtosis of a Gaussian distribution is equal to 3. 
# This indicates that in comparison with the Gaussian distribution, these distributions
# have "thicker tails" - i.e. more probability mass in the tails. 


# Plotting the series
par(mfrow=c(3,1))
plot.ts(m_AA, main = "Return - AA")
plot.ts(m_AXP, main = "Return - AXP")
plot.ts(m_BA, main = "Return = BA")
# It seems that for all return series heteroskedasticity seem to be present. 
# It seems that volatility of returns is not constant over time.
# If we were in the opposite case to heteroskedasticity we would expect 
# the variance to be constant over time. This do however not seem to be the case.
# The graphs shows volatility of the return series in somehow is persistent. 
# If there is a period of high volatility this seems to effect the volatility
# of the upcoming periods. 
# Further it should be noticed that the volatility across the stocks seems to be highly correlated. 




#############################################
##########      Question 2      #############
#############################################

# In this exercise I need to estimate the model of Tse and Tsui (2002) on the three series of returns. 
# I though have a problem with this exercise due to the fact that my previously function
# only can estimate the model for p=2. So therefore I will rewrite the function, so that the function
# computes Psi by vector/matrix notation. 
# I have chosen to rewrite the code for estimation in the following section. This code can take 
# several series of asset return. I have though observed that it seems slower with regards
# to estimating the model. This could potentially be due to the fact that in the computation
# of Psi we invert a matrix I have chosen to call B (see code below) in a loop. These inversions 
# takes a lot of computational power. 

# Note that the following functions follows a similar procedure as from problem set 6. 
tse_tsui_Filter <- function(m_Eta, d_a, d_b, m_Q, d_M) {

    # Main purpose of this is to compute correlation component of log-likelihood 
    # and correlation matrices for each time period t. 
    # The function is used in estimation function in order to estimate Tse and
    # Tsui (2002) model. 
    # Inputs: 
        # m_Eta: Matrix with elements: y_j,t / sigma_j,t for j=1,...,p
        # d_a and d_b: Present in updating equation for R (correlation matrix)
        # m_Q: used for initialization of R
        # d_M: number of periods, where R is constant.

    i_P = ncol(m_Eta)
    i_T = nrow(m_Eta)

    # initialize storage arrays and matrices
    a_R = array(0, dim = c(i_P, i_P, i_T))
    a_D <- array(0, dim = c(i_P, i_P, i_T))
    a_Psi <- array(0, dim = c(i_P, i_P, i_T))

    ## initialization at the unconditional cor
    m_R <- m_Q
    a_R[,, 1] = m_R

    # Initialize LLK
    d_LLK <- t(m_Eta[1, ]) %*% solve(a_R[,, 1]) %*% m_Eta[1, ] - t(m_Eta[1, ]) %*% m_Eta[1, ] + log(det(a_R[,, 1]))
    

    for (i_t in 2:d_M) {
        
        # R matrix, which is constant up to and including period M
        a_R[,,i_t] <- m_R

        # augment the likelihood value (correlation component, slide 35, lecture 10)
        d_LLK <- d_LLK + t(m_Eta[i_t, ]) %*% solve(a_R[,, i_t]) %*% m_Eta[i_t, ] - t(m_Eta[i_t, ]) %*% m_Eta[i_t, ] + log(det(a_R[,, i_t]))
    }


   # Loop for time M up till period T, where R_t is time varying
    for (i_t in seq(i_M + 1, i_T)) {

        # Building the simulation of Psi by matrices and vectors. 
        # This follows the procedure of Tse and Tsui (2002)
        m_B <- matrix(data = 0, nrow = i_P, ncol = i_P, byrow = TRUE)
        diag(m_B) <- (colSums(m_Eta[seq(i_t - 1 ,i_t - i_M), ]^2))
        m_B <- sqrtm(m_B)
        m_E <- t(m_Eta[seq(i_t - 1 ,i_t - i_M) , ]) %*% m_Eta[seq(i_t - 1 ,i_t - i_M) , ]
        a_Psi[,, i_t - 1] <- solve(m_B) %*% m_E %*% solve(m_B)

        # Updating R
        a_R[,,i_t] <- (1 - d_a - d_b) * m_R + d_a * a_Psi[,, i_t - 1] + d_b * a_R[,,i_t - 1]

        # augment the likelihood value (correlation component, slide 35, lecture 10)
        d_LLK <- d_LLK + t(m_Eta[i_t, ]) %*% solve(a_R[,, i_t]) %*% m_Eta[i_t, ] - t(m_Eta[i_t, ]) %*% m_Eta[i_t, ] + log(det(a_R[,, i_t]))
    }

    # Output 
    l_Out = list()
    l_Out[["d_LLK"]] = -0.5 * d_LLK
    l_Out[["a_R"]] = a_R

    return(l_Out)

}


# Objective function which we would like to optimize 
objective <- function(m_Eta, v_Par, m_Q, d_M) {

    Filter = tse_tsui_Filter(m_Eta, v_Par[1], v_Par[2], m_Q, d_M)
    d_NLLK = -as.numeric(Filter$d_LLK)
    return(d_NLLK)

}


Estimate_tse_tsui <- function(m_Y, d_M) {

    # The main purpose of this function is to estimate model of Tse and Tsui (2002). 
    # Input: 
        # m_Y: Series of return
        # d_M: We have that R is constant up to and including time-period M
    # It outputs a list with several results. 
        # d_LLK: Log-likelihood
        # m_Coef: Estimated alpha and beta
        # v_Par: Estimated a and b, which is present in the updating equation for R. 
        # m_Sigma: Volatilities
        # a_R: Correlation matrices for each t
        # m_Eta: Residuals
        # BIC: Bayesian information criteria

    ## estimate the marginal GARCH models
    require(rugarch) #loading the rugarch package
    require(Rsolnp)  #loading the Rsolnp package

    #Marginal garch specifications
    SpecGARCH = ugarchspec(mean.model = list(armaOrder = c(0, 0)))
    m_spec <- multispec( replicate(ncol(m_Y), SpecGARCH) )

    # Fitting garch models
    l_Fit_univariate <- multifit(multispec = m_spec, data = m_Y)

    # Computing residuals
    m_Eta <- residuals(l_Fit_univariate, standardize = TRUE)

    #initial parameters
    v_Par = c(0.003, 0.995)

    # Setting empirical correlation from data
    m_Q <- cor(m_Y)

    # Maximizing likelihood
    optimizer <- solnp(v_Par, objective, 
                        ineqfun = function(v_Par, ...) {
                        sum(v_Par)
                        },
                        ineqLB = 1e-4, ineqUB = 0.999,
                        LB = c(1e-4, 1e-4), UB = c(0.999, 0.999),
                        m_Eta = m_Eta, m_Q = m_Q, d_M = d_M)

    # Extract parameters
    v_Par = optimizer$pars

    #Filter the dynamic correlation using the estimated parameters
    Filter = tse_tsui_Filter(m_Eta, v_Par[1], v_Par[2], m_Q, d_M)

    #extract univariate volatilities
    m_Sigma = sigma(l_Fit_univariate)

    #extract univariate estimated parameters
    m_Coef = coef(l_Fit_univariate)

    #compute the likelihood of the volatility  part
    d_LLK_V = as.numeric(likelihood(l_Fit_univariate))
    d_LLK_V <- sum(d_LLK_V)

    # Extracting likelihood correlation part
    LLK_cor <- -tail(optimizer$values, 1)

    # Computing total likelihood
    d_LLK <- LLK_cor + d_LLK_V
    
    # Correlation matrices
    a_R = Filter[["a_R"]]
   
    # BIC
    BIC = log(i_T) * 8 - 2 * d_LLK

    l_Out = list()

    #output the results
    l_Out[["d_LLK"]] = d_LLK
    l_Out[["m_Coef"]] = m_Coef
    l_Out[["v_Par"]] = v_Par
    l_Out[["m_Sigma"]] = m_Sigma
    l_Out[["a_R"]] = a_R
    l_Out[["m_Eta"]] = m_Eta
    l_Out[["BIC"]] = BIC

    return(l_Out)

}

# Estimating the model of Tse and Tsui (2002).
model_est_5 <- Estimate_tse_tsui(m_Y = m_returns, d_M = 5)
model_est_4 <- Estimate_tse_tsui(m_Y = m_returns, d_M = 4)
model_est_3 <- Estimate_tse_tsui(m_Y = m_returns, d_M = 3)
model_est_2 <- Estimate_tse_tsui(m_Y = m_returns, d_M = 2)
model_est_1 <- Estimate_tse_tsui(m_Y = m_returns, d_M = 1)

# Likelihoods of the different estimations, where M varies from 1 to 5. 
llk_5 <- model_est_5$d_LLK
llk_4 <- model_est_4$d_LLK
llk_3 <- model_est_3$d_LLK
llk_2 <- model_est_2$d_LLK
llk_1 <- model_est_1$d_LLK

m_likelihoods = matrix(c(
    llk_5 <- model_est_5$d_LLK,
    llk_4 <- model_est_4$d_LLK,
    llk_3 <- model_est_3$d_LLK,
    llk_2 <- model_est_2$d_LLK,
    llk_1 <- model_est_1$d_LLK), 
    ncol=5, nrow=1 
)
colnames(m_likelihoods) = c("m=5", "m=4", "m=3", "m=2", "m=1")
m_likelihoods
# We observe that we would choose m = 4, 
# because we here observe the highest log-likelihood. 
#            m=5       m=4       m=3       m=2       m=1
# [1,] -3080.884 -3078.818 -3078.964 -3079.288 -3078.892

# Estimates of a and b
estimates_a_b <- model_est_4$v_Par
# We get: 
    # a = 0.08308817
    # b = 0.79731054

# Further model estimates
estimates_4 <- model_est_4$m_Coef
#               [,1]         [,2]         [,3]
#mu     2.162549e-16 3.722759e-15 2.963473e-15
#omega  7.085599e-03 2.728456e-01 1.170451e-02
#alpha1 8.613071e-02 2.536902e-01 1.802499e-02
#beta1  9.023196e-01 7.415826e-01 9.806518e-01



#############################################
##########      Question 3      #############
#############################################
# Fitting a DCC model using the rmgarch. 
# For this purpose I use the "dccfit" command. 
# First specifying the DCC-GARCH specification prior to fitting
# by using the "dccspec" command. 

# Univariate GARCH specifications
SpecGARCH = ugarchspec(mean.model = list(armaOrder = c(0, 0))) 
# Univariate multiple GARCH specification
m_spec <- multispec( replicate(ncol(m_returns), SpecGARCH) )
# DCC-GARCH specification
dcc_spec <- dccspec(m_spec, model = c("DCC"), distribution = c("mvnorm"))
# Estimating DCC model
dcc_fit <- dccfit(dcc_spec, data = m_returns)
# Obtaining following estimates for a and b: 
#              Estimate  Std. Error   t value Pr(>|t|)
#[Joint]dcca1  0.027362    0.012926  2.116857 0.034272
#[Joint]dccb1  0.946000    0.017155 55.143537 0.000000

# Estimated log likelihood:
dcc_likelihood <- likelihood(dcc_fit)
# Log-Likelihood       :  -3078.857

# Comparing by collecting log-likelihoods into matrix
comparison_likelihood <- matrix(c(
    dcc_likelihood,
    tse_tsui_likelihood = model_est_4$d_LLK),
    ncol=2, nrow=1 
)
colnames(comparison_likelihood) = c("dcc LLK", "Tse Tsui LLK")
# I observe that the estimated likelihoods of the estimations are very close. 
# I though observe that the Tse Tsui estimation has the highest likelihood. 
#  comparison_likelihood
#        dcc LLK   Tse Tsui LLK
# [1,] -3078.857      -3078.818



#############################################
##########      Question 4      #############
#############################################
# Comparing the estimated the correlations of the two models; DCC / Tse and Tsui. 

# First I access the correlations from the DCC estimation
ddc_cor <- rcor(dcc_fit)
# dcc_cor is an array with correlation matrices. 
# The correlation between the different stocks are on the off-diagonal. 
dcc_AA_AXP <- ddc_cor[2, 1, ]
dcc_AA_BA <- ddc_cor[3, 1, ]
dcc_AXP_BA <- ddc_cor[3, 2, ]

# Correlation from Tse and Tsui model
tse_cor <- model_est_4$a_R
tse_AA_AXP <- tse_cor[2, 1, ]
tse_AA_BA <- tse_cor[3, 1, ]
tse_AXP_BA <- tse_cor[3, 2, ]

# Plotting the different correlations of the two models
par(mfrow=c(2,1))
plot.ts(dcc_AA_AXP, main = "DCC model")
lines(dcc_AA_BA, col = "red")
lines(dcc_AXP_BA, col = "purple")
# Adding legend
legend("topleft", legend=c("AA_AXP", "AA_BA", "AXP_BA"),
       col=c("black", "red", "purple"), lty=1:3, cex=0.9)

plot.ts(tse_AA_AXP, main = "Tse and Tsui model")
lines(tse_AA_BA, col = "red")
lines(tse_AXP_BA, col = "purple")
legend("topleft", legend=c("AA_AXP", "AA_BA", "AXP_BA"),
       col=c("black", "red", "purple"), lty=1:3, cex=0.9)


# Next I compare the correlations by plotting the correlations from the 
# different models along with each other. 
par(mfrow=c(3,1))
plot.ts(dcc_AA_AXP, main = "AA AXP")
lines(tse_AA_AXP, col = "red")
legend("topleft", legend=c("DCC", "Tse/Tsui"),
       col=c("black", "red"), lty=1:3, cex=0.9)

plot.ts(dcc_AA_BA, main = "AA BA")
lines(tse_AA_BA, col = "red")
legend("topleft", legend=c("DCC", "Tse/Tsui"),
       col=c("black", "red"), lty=1:3, cex=0.9)

plot.ts(dcc_AXP_BA, main = "AXP BA")
lines(tse_AXP_BA, col = "red")
legend("topleft", legend=c("DCC", "Tse/Tsui"),
       col=c("black", "red"), lty=1:3, cex=0.9)


# Previously I have only been using data for 500 periods.
# Estimating the models considering a time-horizon with 3000 periods.
# I am aware that this might could influence the choice of M. But I  
# chosen to go with M=4 in the following estimation
data("dji30ret")
colnames(dji30ret) # printing tickers
# Defining matrix with data
# I choose the first 3 series from dji30ret and multiply by 100
m_AA <- dji30ret[1:3000, c("AA")] * 100
m_AXP <- dji30ret[1:3000, c("AXP")] * 100
m_BA <- dji30ret[1:3000, c("BA")] * 100

# Removing mean from series
m_AA <- m_AA - mean(m_AA)
m_AXP <- m_AXP - mean(m_AXP)
m_BA <- m_BA - mean(m_BA)

# Combine returns into matrix
m_returns_long <- matrix(data = c(m_AA, m_AXP, m_BA), nrow = length(m_AA), ncol = 3, byrow = FALSE)
colnames(m_returns_long) = c("AA", "AXP", "BA") # setting names for colunms

# Tse and Tsui estimation with m=4 (this might take a while with 3000 observations)
model_est_long <- Estimate_tse_tsui(m_Y = m_returns_long, d_M = 4)
# Correlation from Tse and Tsui model
tse_cor_long <- model_est_long$a_R
tse_AA_AXP <- tse_cor_long[2, 1, ]
tse_AA_BA <- tse_cor_long[3, 1, ]
tse_AXP_BA <- tse_cor_long[3, 2, ]

# DCC estimation
dcc_fit_long <- dccfit(dcc_spec, data = m_returns_long)
# Correlation from DCC model
ddc_cor_long <- rcor(dcc_fit_long)
dcc_AA_AXP <- ddc_cor_long[2, 1, ]
dcc_AA_BA <- ddc_cor_long[3, 1, ]
dcc_AXP_BA <- ddc_cor_long[3, 2, ]

# I compute the estimated volatility of the different stocks, 
# because I want to include it the grid of plots with correlations. 
sigma_est <- model_est_long$m_Sigma
AA_sigma <- sigma_est[, 1]
AXP_sigma <- sigma_est[, 2]
BA_sigma <- sigma_est[, 3]

# Next I compare the correlations by plotting the correlations from the 
# different models along with each other. 
par(mfrow=c(4,1))
plot.ts(dcc_AA_AXP, main = "AA AXP correlation")
lines(tse_AA_AXP, col = "red")
legend("topleft", legend=c("DCC", "Tse/Tsui"),
       col=c("black", "red"), lty=1:3, cex=0.9)

plot.ts(dcc_AA_BA, main = "AA BA correlation")
lines(tse_AA_BA, col = "red")
legend("topleft", legend=c("DCC", "Tse/Tsui"),
       col=c("black", "red"), lty=1:3, cex=0.9)

plot.ts(dcc_AXP_BA, main = "AXP BA correlation")
lines(tse_AXP_BA, col = "red")
legend("topleft", legend=c("DCC", "Tse/Tsui"),
       col=c("black", "red"), lty=1:3, cex=0.9)

plot.ts(AA_sigma, main = "Volatilities")
lines(AXP_sigma, col = "red")
lines(BA_sigma, col = "purple")
legend("topleft", legend=c("AA", "AXP", "BA"),
       col=c("black", "red", "purple"), lty=1:3, cex=0.9)




#############################################
##########      Question 5      #############
#############################################
# Conclusions about the estimated volatilities and correlations. 

# If we start by considering the estimated correlations for T=3000, we have
# that the estimated correlations from the Tse and Tsui model seem to be more
# stable around its mean compared to the estimated
# correlation of the DCC model. 
# Especially in the periods where the DCC model estimates high or low values for
# the correlation (compared to its mean), we see that the estimated correlation
# of the Tse and Tsui model seem to keep fluctuating around its mean.
# To further emphasis this, it is possible to compare the std. deviation of the
# correlations. 
# Below we observe that the std. deviation of the estimated correlations of the
# Tse and Tsui model are smaller than of the std. deviation the estimated 
# correlations of the DCC model. 
# Comparison of standard for correlation of AA/AXP
comparison_sd_cor_AA_AXP <- matrix(c(
    sd_tse_AA_AXP <- sd(tse_AA_AXP),
    sd_dcc_AA_AXP = sd(dcc_AA_AXP)),
    ncol=2, nrow=1 
)
colnames(comparison_sd_cor_AA_AXP) = c("TSE/TSUI", "DCC")
#       TSE/TSUI       DCC
#[1,] 0.06741146 0.1173605

# When considering the plot from question 4 (empirical part) where volatilies(sigma)
# of the different assets are illustrated in the bottom panel, it can be observed that 
# in periods with higher volatility the correlation between assets seems
# to react by increasing correlation. 
# It should though be noticed that this relation between correlation and volatility 
# seems to be mostly present for the DCC estimated correlations. 






