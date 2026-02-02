
# Load required libraries
library(ctsmTMB)

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

### Create NN-embedded ctsmTMB model

# Neural network settings
nn_settings <- list(
    Na = list(
        input = c("ca"),
        no_neurons = NA,
        act_fun = NA
    ),
    Nd = list(
        input = c("ca", "cd", "F", "Q"),
        no_neurons = c(8),
        act_fun = "tanh"
    )
)

x0 <- c(2.0, 0.5, 0.5) # Initial state
P0 <- diag(c(0.1, 0.1, 0.1)) # Initial state covariance
input_names <- c("F", "c_gina", "F_gina", "Q")   # Input variables

sofun()
model <- create_ctsm_model(x0, P0, 
                                input_names,    
                                incl_intermediate_storage = TRUE,
                                nn_settings)

#### 


# Create random data for testing
df <- data.frame(
    t = seq(0, 10, by = 0.1),
    F = runif(101, 0.5, 1.5),
    c_gina = runif(101, 1.0, 3.0),
    F_gina = runif(101, 0.5, 1.5),
    Q = runif(101, 0.1, 0.5),
    yga = runif(101, 20, 80)
)

# Get the likelihood function

nll_fn <- model$likelihood(df, 
            method = "laplace",
            ode.solver = "euler", 
            ode.stepsize = 0.01
)

par <- model$getParameters("type" = "free")[,"initial"]
par

# Evaluate the negative log-likelihood
nll_value <- nll_fn$fn(par)
nl_grad <- nll_fn$gr(par)
nll_value
nl_grad
