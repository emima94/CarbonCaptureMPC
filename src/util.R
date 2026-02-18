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
predict_and_plot <- function(model, df, par, x0, P0, ode_solver = "rk4", ode_stepsize = 10/3600) {

    df_pred <- df
    #df_pred$yga <- NA  # Set yga to NA to indicate it's not observed during prediction


    pred <- model$predict(
        data = df_pred,
        method = "ekf",
        ode_solver = ode_solver,
        ode_stepsize = ode_stepsize,
        pars = par,
        initial.state = list(x0, P0)
        #estimate.initial.state = TRUE
    )

    pred <- post_process_pred(model, df, pred, par_free = par)

    t_pred <- pred$observations$t.j
    y_pred <- pred$observations$yga

    ### PLOT ###

    par(mfrow = c(5,1))
    # yga plot, with CI, with data points
    plot(df$t, df$yga, 
        type = "p", col = "red", xlab = "Time (h)", ylab = "yga (%)", 
        main = "Initial prediction", ylim = c(0, 4))
    lines(t_pred, y_pred, col = "blue")
    lines(t_pred, pred$observations$yga_CImax, col = "lightblue", lty = 2)
    lines(t_pred, pred$observations$yga_CImin, col = "lightblue", lty = 2)
    legend("topright", legend = c("Data", "Prediction"), col = c("red", "blue"), lty = c(NA, 1), pch = c(1, NA))
    # ca and cd plots, predictions only, with CI
    plot(t_pred, pred$states$ca, type = "l", col = "blue", 
        xlab = "Time (h)", ylab = "ca", 
        main = "ca prediction",
        ylim = c(0,3))
    lines(t_pred, pred$states$ca_CImax, col = "lightblue", lty = 2)
    lines(t_pred, pred$states$ca_CImin, col = "lightblue", lty = 2)
    lines(t_pred, pred$states$cd, type = "l", col = "red", xlab = "Time (h)", ylab = "cd", main = "cd prediction", ylim = c(0,3))
    lines(t_pred, pred$states$cd_CImax, col = "pink", lty = 2)
    lines(t_pred, pred$states$cd_CImin, col = "pink", lty = 2)
    # Na and Nd plots, include data for Na
    plot(df$t, df$Na, type = "p", col = "black", ylim = c(300, 1000))
    lines(t_pred, pred$algebraics$Na, col = "blue", lty = 2)
    lines(t_pred, pred$algebraics$Nd, col = "red", lty = 2)
    legend("topright", legend = c("Na data", "Na pred", "Nd pred"), col = c("black", "blue", "red"), lty = c(NA, 2, 2), pch = c(1, NA, NA))
    
    plot(df$t, df$F, type = "l", col = "green", ylim = c(0.2, 0.4))
    
    plot(df$t, df$Q, type = "l", col = "purple", ylim = c(22, 34))
    
    return(pred)
}