compute_ss_value <- function(df, var, k = 10) {
    # df: data frame with time series data
    # var: variable name (string) for which to compute steady state
    # k: number of time points to average over for steady state
    start_ss <- mean(df[[var]][1:k])
    end_ss <- mean(tail(df[[var]], k))
    return(c(start_ss, end_ss))
}

get_ss_data <- function(df, train_idxs, vars_ss, k = 10) {
    ss_list <- list()
    for (var in vars_ss) {
        ss_values <- c()
        for (i in train_idxs) {
            ss_i <- compute_ss_value(df[[i]], var, k)
            ss_values <- c(ss_values, ss_i)
        }
        ss_list[[var]] <- ss_values
    }
    df_ss <- as.data.frame(ss_list)
    return(df_ss)
}

# Test the initial parameter guess with predictions on dataset i
predict_and_plot <- function(model, df, i, par, x0, P0, ode_solver = "rk4", ode_stepsize = 10/3600) {

    pred <- model$predict(
        data = df[[i]],
        method = "ekf",
        ode_solver = ode_solver,
        ode_stepsize = ode_stepsize,
        pars = par,
        initial.state = list(x0[[i]], P0)
    )

    t_pred <- pred$observations$t.j
    y_pred <- pred$observations$yga

    plot(df[[i]]$t, df[[i]]$yga, 
        type = "p", col = "red", xlab = "Time (h)", ylab = "yga (%)", 
        main = sprintf("Data series %d: Initial prediction", i),
        ylim = c(0,15))
    lines(t_pred, y_pred, col = "blue")
    legend("topright", legend = c("Data", "Prediction"), col = c("red", "blue"), lty = c(NA, 1), pch = c(1, NA))
}