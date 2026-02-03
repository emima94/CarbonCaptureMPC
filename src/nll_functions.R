# Aggregate NLL function over all training series
nll_fn <- function(par, nll_fun_list) {
    nll_total <- 0
    for (nll_fn_i in nll_fun_list) {
        nll_total <- nll_total + nll_fn_i$fn(par)
    }
    return(nll_total)
}
nll_gr <- function(par, nll_fun_list) {
    nll_grad_total <- rep(0, length(par))
    for (nll_fn_i in nll_fun_list) {
        nll_grad_total <- nll_grad_total + nll_fn_i$gr(par)
    }
    return(nll_grad_total)
}
