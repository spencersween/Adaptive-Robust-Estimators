rm(list = ls())
set.seed(123)

# Load Packages Here
if (!require(numDeriv)) {
  install.packages("numDeriv")
}
if (!require(tidyverse)) {
  install.packages("tidyverse")
}
if (!require(RColorBrewer)) {
  install.packages("tidyverse")
}
library(numDeriv)
library(tidyverse)

# Define Functions Here
tau_finder = function(X, Y, b) {

    r = Y - X %*% b
    n = dim(X)[1]
    d = dim(X)[2]
    
    tau_function = function(x) {
        (sum(pmin(r ^ 4, x)) / (x * (d + log(n)))) - 1
    }

    tau_min = min(r ^ 4)
    tau_max = sum(r ^ 4)

    root = try(uniroot(tau_function, c(tau_min, tau_max))$root, silent = TRUE)

    if (class(root) == "try-error") {
        tau = (sum(r ^ 4) / (d + log(n))) ^ (1 / 4)
    }
    
    else {
        tau = (root) ^ (1 / 4)
    }

    return(tau)

}

huber_loss_function = function(b, X, Y, index, tau = 1.345) {
    r = Y[index] - X[index, ] %*% b
    true_false = (abs(r) <= tau)
    loss = ifelse(true_false, ((r ^ 2) / 2), ((tau * (abs(r)) - (tau / 2))))
    return(sum(loss))
}

capped_loss_function = function(b, X, Y, index, tau = 1.345) {
    r = Y[index] - X[index, ] %*% b
    true_false = (abs(r) <= tau)
    loss = ifelse(true_false, ((r ^ 2) / 2), 0)
    return(sum(loss))
}

create_plot = function(data, variable) {

    # Filter data for the specified variable
    data_filtered = data %>% dplyr::filter(Variable == variable)
    
    # Calculate the average of the estimates
    avg_estimate = mean(data_filtered$Estimate)
    
    # Create the plot
    Plot = ggplot2::ggplot(data_filtered, ggplot2::aes(x = Iteration, y = Estimate)) +
        ggplot2::geom_line(size = 1.5, color = "#2C3E50") +
        ggplot2::geom_hline(yintercept = avg_estimate, linetype = "dashed", color = "black", size = 2) +
        ggplot2::labs(
        title = paste("Iteration Path of", variable, "estimates"),
        x = "Iteration",
        y = "Estimate",
        caption = paste("Black dashed line represents the average of", variable, "estimates")
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggthemes::theme_tufte() +
        ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
        axis.title = ggplot2::element_text(size = 16),
        legend.position = "none"
        )

    return(Plot)
}

so_robust_regression = function(X, Y, type = "Huber", iter = 1000, lr = 0.1000, m_factor = 1, print_var) {

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
        tau_b[i] = tau_i
        use_loss = ifelse(type == "Huber", huber_loss_function, capped_loss_function)
        G = jacobian(use_loss, beta_i, X = X, Y = Y, index = index_i, tau = tau_i)
        H = hessian(use_loss, beta_i, X = X, Y = Y, index = index_i, tau = tau_i)
        beta_i = beta_i - lr * solve(H) %*% t(G)
        if(i > burn) {
            beta_b[i-burn,] = beta_i
        }
    }

    # Create Plot
    beta_b_df = as.data.frame(beta_b)
    names(beta_b_df) = colnames(X)
    beta_b_df$Iteration <- (burn+1):(burn+iter)
    beta_b_long <- tidyr::pivot_longer(beta_b_df, -Iteration, names_to = "Variable", values_to = "Estimate")
    Plot = create_plot(beta_b_long, print_var)

    # Create Table
    estimates = apply(beta_b, 2, mean)
    ses = apply(beta_b, 2, sd) * se_adjust
    tstat = estimates / ses
    cil = estimates - 1.96 * ses
    ciu = estimates + 1.96 * ses
    reject = ifelse(cil > 0 | ciu < 0, 1, 0)
    Table = tibble(
        `Variable` = colnames(X),
        `Loss` = ifelse(type == "Huber", "AHR", "ACLS"),
        `Estimate` = sprintf("%.2f", estimates), 
        `SE` = sprintf("%.2f", ses), 
        `CI_Lower` = sprintf("%.2f", cil), 
        `CI_Upper` = sprintf("%.2f", ciu), 
    )
    return(list(Table = Table, Plot = Plot))
}

data = read_csv("/Users/spencersween/Downloads/blimpo.csv") %>% mutate(Treatment = d)
X = data %>% dplyr::select(Treatment) %>% mutate(Intercept = rep(1, 789)) %>% as.matrix()
Y = data %>% dplyr::select(y) %>% as.matrix()
Results = so_robust_regression(X = X, Y = Y, type = "Huber", iter = 10000, lr = 0.0100, m_factor = 1, print_var = "Treatment")
print(Results$Table)
Results$Plot