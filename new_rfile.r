rm(list = ls())
# This file contains the code for the functions used in the robust regression simulation study.

############################################
##### Load Packages Here #####
############################################

if (!require(numDeriv)) {
  install.packages("numDeriv")
}
if (!require(tidyverse)) {
  install.packages("tidyverse")
}

library(numDeriv)
library(tidyverse)


############################################
##### Define Functions Here #####
############################################

# Define: Tuning Parameter Search Function
tau_finder = function(X, Y, b) {
    
    # Calculate residuals
    r = Y - X %*% b
    
    # Get the number of observations (n) and variables (d) from X
    n = dim(X)[1]
    d = dim(X)[2]
    
    # Define a helper function tau_function that takes an input x
    tau_function = function(x) {
        # Return the value of the function for the given x
        (sum(pmin(r ^ 4, x)) / (x * (d + log(n)))) - 1
    }

    # Calculate the minimum and maximum values for tau
    tau_min = min(r ^ 4)
    tau_max = sum(r ^ 4)

    # Try to find the root of the tau_function in the interval [tau_min, tau_max]
    root = try(uniroot(tau_function, c(tau_min, tau_max))$root, silent = TRUE)

    # If finding the root failed (i.e., root is a "try-error")
    if (class(root) == "try-error") {
        # Calculate tau using a fallback formula
        tau = (sum(r ^ 4) / (d + log(n))) ^ (1 / 4)
    }
    
    # If finding the root succeeded
    else {
        # Calculate tau using the root
        tau = (root) ^ (1 / 4)
    }

    # Return the value of tau
    return(tau)

}

# Define: Robust Loss Function
robust_loss_function = function(b, X, Y, index, tau, type = "Huber") {

    # Calculate the residuals
    r = Y[index] - X[index, ] %*% b
    
    # Determine if the residuals are within the threshold tau
    true_false = (abs(r) <= tau)
    
    # Calculate the loss based on the residuals, threshold, and type

    if(type == "Huber") {
        loss = sum(ifelse(true_false, ((r ^ 2) / 2), ((tau * (abs(r)) - (tau / 2)))))
    }
    if(type == "Capped") {
        loss = sum(ifelse(true_false, ((r ^ 2) / 2), 0))
    }
    
    # Return the sum of the loss
    return(loss)
}

jacobian_robust_loss_function = function(b, X, Y, index, tau, type = "Huber") {

    # Calculate the residuals
    r = Y[index] - X[index, ] %*% b
    
    # Determine if the residuals are within the threshold tau
    true_false = (abs(r) <= tau)

    # Calculate the gradient based on the residuals, threshold, and type

    if(type == "Huber") {
        grad = t(X[index,]) %*% ifelse(true_false, -r, tau * sign(-r))
    }
    if(type == "Capped") {
        grad = t(X[index,]) %*% ifelse(true_false, -r, 0)
    }
    
    # Return the gradient
    return(grad)
}

hessian_robust_loss_function = function(b, X, Y, index, tau) {

    # Calculate the residuals
    r = Y[index] - X[index, ] %*% b
    
    # Determine if the residuals are within the threshold tau
    true_false = (abs(r) <= tau)

    # Calculate the hessian
    hessian = t(X[index,]) %*% diag(as.numeric(true_false)) %*% X[index,]
    
    # Return the hessian
    return(hessian)
}


# Define: Stochastic Optimization Estimation and Inference Procedure
so_robust_regression = function(X, Y, type = "Huber", iter = 200, lr = 0.10, m_factor = 1) {

        n = dim(X)[1]
        d = dim(X)[2]

        m = m_factor * n
        burn = 1 + round(log(0.01) / log(1 - lr))
        se_adjust = sqrt(m / n) * sqrt((1 - (1 - lr) ^ 2) / lr ^ 2)

        beta_b = matrix(NA, nrow = iter, ncol = d)
        tau_b = matrix(NA, nrow = iter, ncol = 1)

        beta_i = lm.fit(X,Y)$coefficients
        for(i in 1:(iter+burn)) {
                index_i = sample(1:n, size = n, replace = TRUE)
                tau_i = tau_finder(X[index_i, ], Y[index_i], beta_i)

                g = jacobian_robust_loss_function(b = beta_i, X = X, Y = Y, index = index_i, tau = tau_i, type = type)
                H = hessian_robust_loss_function(b = beta_i, X = X, Y = Y, index = index_i, tau = tau_i)

                # Use Damped Semi-smooth Newton Method
                H_inverse = try(solve(H), silent = TRUE)
                if (class(H_inverse)[1] == "try-error") {
                    Update = lr * g
                }
                else {
                    Update = lr * H_inverse %*% g
                }
                
                beta_i = beta_i - Update
                if(i > burn) {
                        tau_b[i-burn,] = tau_i
                        beta_b[i-burn,] = beta_i
                }
        }

        # Table
        estimates = apply(beta_b, 2, mean)
        ses = apply(beta_b, 2, sd) * se_adjust
        tstat = estimates / ses
        cil = estimates - 1.96 * ses
        ciu = estimates + 1.96 * ses
        reject = ifelse(cil > 0 | ciu < 0, 1, 0)
        Table = tibble(
            `Variable` = colnames(X),
            `Loss` = ifelse(type == "Huber", "AHR", "ACLS"),
            `Tau` = sprintf("%.0f", mean(tau_b)),
            `Estimate` = sprintf("%.2f", estimates), 
            `SE` = sprintf("%.2f", ses), 
            `CI_Lower` = sprintf("%.2f", cil), 
            `CI_Upper` = sprintf("%.2f", ciu), 
        )

        # Return Objects
        estimate = apply(beta_b, 2, mean)
        se = apply(beta_b, 2, sd) * se_adjust
        tau = mean(tau_b)
        return(list(Estimate = estimate, SE = se, Tau = tau, Table = Table))


}



data = read_csv("/Users/spencersween/Downloads/blimpo.csv")
X = data %>% dplyr::select(-y, -d) %>% as.matrix()
D = data %>% dplyr::select(d) %>% as.matrix()
Y = data %>% dplyr::select(y) %>% as.matrix()

Y_tilde = Y - X %*% lm.fit(X, Y)$coefficients
D_tilde = cbind(D - X %*% lm.fit(X, D)$coefficients, rep(1, 789))
colnames(D_tilde) = c("Treatment", "Intercept")
result = so_robust_regression(X = D_tilde, Y = Y_tilde, type = "Huber", iter = 200, lr = 0.30, m_factor = 1)
print(result$Table)

D_2 = cbind(D, rep(1, 789))
colnames(D_2) = c("Treatment", "Intercept")
result = so_robust_regression(X = D_2, Y = Y, type = "Huber", iter = 200, lr = 0.30, m_factor = 1)
print(result$Table)
summary(MASS::rlm(Y~D))
summary(lm(Y~D))

