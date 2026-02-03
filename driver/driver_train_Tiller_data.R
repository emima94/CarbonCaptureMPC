
# Load required libraries
library(ctsmTMB)
library(readxl)   # for reading Excel
library(dplyr)    # for data manipulation
library(readr)    # for reading CSV
library(nleqslv)  # for solving nonlinear equations

# Seed random number generator for reproducibility
set.seed(123)

# Source helper functions
# Source all functions in src/
sofun <- function() {
    src_files <- list.files("src",pattern="*.R",full.names=TRUE)
    for (f in src_files) {
        source(f)
    }
}
sofun()

## ----------------------------------- ##
#### Load data ####
## ----------------------------------- ##

df <- load_data_Tiller()


# Steady state data variables
train_idxs <- 1:10  # Indices of training data series
input_vars_ss <- c("F", "Q", "ca", "cd", "cgina", "Fgina")
df_ss_inputs <- get_ss_data(df, train_idxs, input_vars_ss, k = 10)

output_vars_ss <- c("Na")
df_ss_outputs <- get_ss_data(df, train_idxs, output_vars_ss, k = 10)

## ----------------------------------- ##
#### Model training / MLE ####
## ----------------------------------- ##

# Compute/choose mean and standard deviation of NN inputs and outputs
ca_vec <- unlist(sapply(df, function(d) d$ca))
ca_mean <- mean(ca_vec)
ca_sd <- sd(ca_vec)

cd_vec <- unlist(sapply(df, function(d) d$cd))
cd_mean <- mean(cd_vec)
cd_sd <- sd(cd_vec)

F_vec <- unlist(sapply(df, function(d) d$F))
F_mean <- mean(F_vec)
F_sd <- (max(F_vec) - min(F_vec)) / 2
Q_vec <- unlist(sapply(df, function(d) d$Q))
Q_mean <- mean(Q_vec)
Q_sd <- (max(Q_vec) - min(Q_vec)) / 2
Na_vec <- unlist(sapply(df, function(d) d$Na))
Na_mean <- mean(Na_vec)
Na_sd <- sd(Na_vec)    
cat(sprintf("Chosen ca mean: %.2f, sd: %.2f\n", ca_mean, ca_sd))
cat(sprintf("Chosen cd mean: %.2f, sd: %.2f\n", cd_mean, cd_sd))
cat(sprintf("Computed F mean: %.2f, sd: %.2f\n", F_mean, F_sd))
cat(sprintf("Computed Q mean: %.2f, sd: %.2f\n", Q_mean, Q_sd))
cat(sprintf("Computed Na mean: %.2f, sd: %.2f\n", Na_mean, Na_sd))

plot(cd_vec)

### Model settings, ctsmTMB model ###
# NN settings
nn_settings <- list(
    Na = list(
        input = list(
            names = c("ca"),
            means = c(ca_mean),
            sds = c(ca_sd)
        ),
        output = list(
            name = "Na",
            mean = Na_mean,
            sd = Na_sd
        ),
        no_neurons = NA,
        act_fun = NA
    ),
    Nd = list(
        input = list(
            names = c("ca", "cd", "F", "Q"),
            means = c(ca_mean, cd_mean, F_mean, Q_mean),
            sds = c(ca_sd, cd_sd, F_sd, Q_sd)
        ),
        output = list(
            name = "Nd",
            mean = Na_mean,
            sd = Na_sd
        ),
        no_neurons = c(8),
        act_fun = "tanh"
    )
)

# Initial state and variance
x0 <- list()
k <- 10  # number of points to average for initial state
for (i in seq_along(df)) {
    x0[[i]] <- c(mean(df[[i]]$ca[1:k]), 
                 mean(df[[i]]$cd[1:k]), 
                 mean(df[[i]]$cd[1:k]))
}
P0 <- diag(c(0.1, 0.1, 0.1))
# input names
input_names <- c("F", "cgina", "Fgina", "Q")
# include intermediate storage compartment
incl_intermediate_storage <- TRUE


# Training settings
train_idxs <- c(1,2,3,4,5,6,7,8,9,10)  # Indices of training data series
method <- "laplace"
ode_solver <- "rk4"
ode_stepsize <- 10 / 3600  # 10 seconds

# Initial parameter guess and lower and upper bounds
model <- create_ctsm_model(x0[[1]], P0, incl_intermediate_storage, nn_settings)
par <- model$getParameters("type" = "free")[,"initial"]
lb <- model$getParameters("type" = "free")[,"lower"]
ub <- model$getParameters("type" = "free")[,"upper"]

# Plot initial predictions for training series
for (i in train_idxs) {
    #predict_and_plot(model, df, i, par, ode_solver, ode_stepsize)
}

#
sofun()


ss_sol <- nlminb(start = par,
                objective = g_ss_fn,
                lower = lb,
                upper = ub,
                control = list(eval.max = 500, 
                                iter.max = 500,
                                trace = 10)
                )

par_init_ss <- ss_sol$par


#predict_and_plot(model, df, 2, par_init_ss, ode_solver, ode_stepsize)



# Compile NLL functions for training series
nll_fun_list <- compile_nll(df, train_idxs, x0, P0, 
                            nn_settings,
                            method,
                            ode_solver,
                            ode_stepsize,
                            incl_intermediate_storage)


# Initial nll
nll_initial <- nll_fn(par, nll_fun_list)
nll_initial_ss_est <- nll_fn(par_init_ss, nll_fun_list)
cat(sprintf("Initial NLL: %.4f\n", nll_initial))
cat(sprintf("Initial NLL (ss est): %.4f\n", nll_initial_ss_est))
# Initial gradient norm
nll_gr_initial <- nll_gr(par_init_ss, nll_fun_list)
grad_norm_initial <- sqrt(sum(nll_gr_initial^2))
cat(sprintf("Initial gradient norm: %.4f\n", grad_norm_initial))

# Parameters, steady state solutions

# Optimize parameters to minimize NLL
fit <- nlminb(start = par_init_ss,
                objective = function(par) nll_fn(par, nll_fun_list),
                gradient = function(par) nll_gr(par, nll_fun_list),
                lower = lb,
                upper = ub,
                control = list(eval.max = 1000, 
                                iter.max = 1000,
                                trace = 10)
                )

par_star <- fit$par
# Predict and plot after training
sofun()
i <- 10
predict_and_plot(model, df, i, par_star, x0, P0, ode_solver, ode_stepsize)

