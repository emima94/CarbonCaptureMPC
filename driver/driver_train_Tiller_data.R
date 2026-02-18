
# Load required libraries
library(ctsmTMB)
library(readxl)   # for reading Excel
library(dplyr)    # for data manipulation
library(readr)    # for reading CSV
library(nleqslv)  # for solving nonlinear equations
library(rlang)   # for expression handling
library(jsonlite) 

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
# i <- 7
# plot(df[[i]]$t, df[[i]]$cga, type = "l", col = "red", xlab = "Time (h)", ylab = "yga (%)", main = sprintf("Data series %d: yga", i))

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

cga_vec <- unlist(sapply(df, function(d) d$cga))
cga_mean <- mean(cga_vec)
cga_sd <- sd(cga_vec)

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

plot(unlist(sapply(df, function(d) d$cga)))

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
ode_stepsize <- 5 / 3600  # 

# Initial parameter guess and lower and upper bounds
sofun()
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

i <- 3
predict_and_plot(model, df[[i]], par_init_ss, x0[[i]], P0, ode_solver, ode_stepsize)



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
                control = list(eval.max = 10000, 
                                iter.max = 10000,
                                trace = 10)
                )

par_star <- fit$par
names(par_star) <- rownames(model$getParameters("type" = "free"))

# Write par est to csv
write.csv(par_star, "results/parameter_estimates_Tiller.csv", row.names = TRUE)
# Write fixed par to csv
par_fixed <- model$getParameters("type" = "fixed")[,"initial"]
names(par_fixed) <- rownames(model$getParameters("type" = "fixed"))
write.csv(par_fixed, "results/parameter_fixed_Tiller.csv", row.names = TRUE)
# Write nn_settings to json
write_json(nn_settings, "results/nn_settings_Tiller.json", pretty = TRUE)




# Evaluate NNs for slices of NN input space
F <- seq(0.2, 0.45, length.out = 50)

Nd_F <- sapply(F, function(F_i) {
    df_input <- data.frame(ca = 2.63, cd = 0.5, F = F_i, Q = 28, cgina = NaN, Fgina = NaN)
    evaluate_algebraics_simple(model, df_input, par_free = par_star)$Nd
})
Nd_F
par(mfrow = c(1,1))
plot(F, Nd_F, type = "l", xlab = "F", ylab = "Nd", main = "Estimated Nd vs F")
Q <- seq(22,34, length.out = 50)
Nd_Q <- sapply(Q, function(Q_i) {
    df_input <- data.frame(ca = 2.63, cd = 0.4, F = 0.32, Q = Q_i, cgina = NaN, Fgina = NaN)
    evaluate_algebraics_simple(model, df_input, par_free = par_star)$Nd
})
plot(Q, Nd_Q, type = "l", xlab = "Q", ylab = "Nd", main = "Estimated Nd vs Q")
ca <- seq(2.0, 3.5, length.out = 50)
Nd_ca <- sapply(ca, function(ca_i) {
    df_input <- data.frame(ca = ca_i, cd = 0.4, F = 0.32, Q = 28, cgina = NaN, Fgina = NaN)
    evaluate_algebraics_simple(model, df_input, par_free = par_star)$Nd
})
plot(ca, Nd_ca, type = "l", xlab = "ca", ylab = "Nd", main = "Estimated Nd vs ca")


# Predict and plot after training
sofun()
i <- 8
pred <- predict_and_plot(model, df[[i]], par_star, x0[[i]], P0, ode_solver, ode_stepsize)

# Plot data
cgina_vec <- unlist(sapply(df, function(df) df$cgina))
Fgina_vec <- unlist(sapply(df, function(df) df$Fgina))

plot(cgina_vec, type = "l", col = "red", xlab = "Time (h)", ylab = "cgina", main = sprintf("Data series %d: cgina", i))
plot(Fgina_vec, type = "l", col = "blue", xlab = "Time (h)", ylab = "Fgina", main = sprintf("Data series %d: Fgina", i))

# Simulate with new inputs and trained parameters
#case_name <- "case1"
t0_sim <- 0
tf_sim <- 5
dt_sim <- 30 / 3600  # 30 second intervals
nt <- floor((tf_sim - t0_sim) / dt_sim) + 1

# Compute ygina
Pa <- 101325
R <- 8.314
Ta <- 313.15
cgina_tot <- Pa /(R * Ta)

# Values for step changes
cgina_values <- c(4.4, 4.4, 4.2, 4.2, 4.8, 4.8)
Fgina_values <- c(180, 160, 160, 200, 200, 180)
F_values <- c(0.3, 0.2, 0.1, 0.01, 0.5, 0.6)#c(0.3, 0.4, 0.4, 0.22, 0.22, 0.3)
Q_values <- c(28, 28, 24, 24, 31, 31)
n_steps <- length(cgina_values)
step_length <- ceiling(nt / n_steps)
df_list <- list(
    df_case1 <- data.frame(
    t = seq(t0_sim, tf_sim, by = dt_sim),  # 30 second intervals for 5 hours
    yga = NaN,  # NA for yga since we're simulating
    F = 0.3,  # Constant F
    cgina = 4.8,  # Constant cgina
    Fgina = 160,  # Constant Fgina
    Q = 28  # Constant Q
    ),
    df_case2 <- data.frame(
        t = seq(t0_sim, tf_sim, by = dt_sim),  # 30 second intervals for 5 hours
        yga = NaN,  # NA for yga since we're simulating
        F = 0.3,  # Constant F
        cgina = rep(cgina_values, each = step_length)[1:nt],  # Step changes in cgina
        Fgina = rep(Fgina_values, each = step_length)[1:nt],  # Step changes in Fgina
        Q = 28  # Constant Q    
    ),
    df_case3 <- data.frame(
        t = seq(t0_sim, tf_sim, by = dt_sim),
        yga = NaN,  # NA for yga since we're simulating
        F = rep(F_values, each = step_length)[1:nt],  # Step changes in F
        cgina = 4.8,
        Fgina = 160,
        Q = rep(Q_values, each = step_length)[1:nt]  # Step changes in Q
    )

)

for (i in seq_along(df_list)) {
    df_i <- df_list[[i]]
    case_name <- paste0("case", i)

    pred_case_i <- model$predict(
        data = df_i,
        method = "ekf",
        ode_solver = ode_solver,
        ode_stepsize = ode_stepsize,
        pars = par_star,
        initial.state = list(x0[[1]], P0)
        )

    pred_case_i <- post_process_pred(model, df_i, pred_case_i, par_free = par_star)

    ylims <- list(
        yga = c(-10, 16),
        ca = c(-1, 3.5),
        #cd = c(0, 2.8),
        Na = c(400, 900),
        Nd = c(400, 900),
        F = c(0.0, 0.5),
        cgina = c(4.0, 5.0),
        Fgina = c(150, 210),
        Q = c(20, 35)
    )

    ## Outputs ##
    scale <- 7
    ratio <- 1.0
    pdf(sprintf("figures/simulations/%s_outputs.pdf", case_name), width = scale, height = ratio * scale)
    par(mfrow = c(3,1))
    # yga simulation, incl. CI (polygon), incl ygina input as reference (dashed line)
    plot(pred_case_i$observations$t.j, pred_case_i$observations$yga, type = "n", xlab = "Time (h)", ylab = "yga (%)", main = "Simulated yga with trained model", ylim = ylims$yga)
    polygon(c(pred_case_i$observations$t.j, rev(pred_case_i$observations$t.j)),
            c(pred_case_i$observations$yga_CImin, rev(pred_case_i$observations$yga_CImax)),
            col = rgb(0, 0, 1, alpha = 0.2), border = NA)
    lines(pred_case_i$observations$t.j, pred_case_i$observations$yga, type = "l", col = "blue")
    lines(df_i$t, df_i$cgina / cgina_tot * 100, type = "l", col = "red", lty = 2)
    legend("topright", legend = c("yga", "ygina (input)"), col = c("blue", "red"), lty = c(1, 2)) 
    # simulated states, incl. CI (polygon)

    plot(pred_case_i$states$t.j, pred_case_i$states$ca, type = "n", xlab = "Time (h)", ylab = "ca", main = "Simulated ca", ylim = ylims$ca)
    polygon(c(pred_case_i$states$t.j, rev(pred_case_i$states$t.j)),
            c(pred_case_i$states$ca_CImin, rev(pred_case_i$states$ca_CImax)),
            col = rgb(0, 0, 1, alpha = 0.2), border = NA)
    polygon(c(pred_case_i$states$t.j, rev(pred_case_i$states$t.j)),
            c(pred_case_i$states$cd_CImin, rev(pred_case_i$states$cd_CImax)),
            col = rgb(1, 0, 0, alpha = 0.2), border = NA)
    #lines(pred_case_i$states$t.j, pred_case_i$states$cd_CImin, type = "l", col = rgb(1, 0, 0, alpha = 0.5), lty = 2)
    #lines(pred_case_i$states$t.j, pred_case_i$states$cd_CImax, type = "l", col = rgb(1, 0, 0, alpha = 0.5), lty = 2)
    lines(pred_case_i$states$t.j, pred_case_i$states$ca, type = "l", col = "blue")
    lines(pred_case_i$states$t.j, pred_case_i$states$cd, type = "l", col = "red")
    legend("topright", legend = c("ca", "cd"), col = c("blue", "red"), lty = 1)
    # simulated Na and Nd
    plot(pred_case_i$states$t.j, pred_case_i$algebraics$Na, type = "l", col = "blue", xlab = "Time (h)", ylab = "Na", main = "Simulated Na", ylim = ylims$Na)
    lines(pred_case_i$states$t.j, pred_case_i$algebraics$Nd, type = "l", col = "red", xlab = "Time (h)", ylab = "Nd", main = "Simulated Nd", ylim = ylims$Nd)
    legend("topright", legend = c("Na", "Nd"), col = c("blue", "red"), lty = 1)
    dev.off()
    ## Inputs ##
    pdf(sprintf("figures/simulations/%s_inputs.pdf", case_name), width = scale, height = ratio * scale)
    par(mfrow = c(4,1))
    plot(df_i$t, df_i$cgina, type = "l", col = "red", xlab = "Time (h)", ylab = "cgina", main = "Input cgina", ylim = ylims$cgina)
    plot(df_i$t, df_i$Fgina, type = "l", col = "green", xlab = "Time (h)", ylab = "Fgina", main = "Input Fgina", ylim = ylims$Fgina)
    plot(df_i$t, df_i$F, type = "l", col = "blue", xlab = "Time (h)", ylab = "F", main = "Input F", ylim = ylims$F)
    plot(df_i$t, df_i$Q, type = "l", col = "purple", xlab = "Time (h)", ylab = "Q", main = "Input Q", ylim = ylims$Q)
    dev.off()
}

# Open-loop MPC simulation with trained model.


