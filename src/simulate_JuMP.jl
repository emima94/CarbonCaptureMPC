function euler_simulation(ffun, hfun, x0, U, D, p, t_vec, dt_sim)
    nx = length(x0)
    X = zeros(nx, length(t_vec))
    X[:,1] = x0
    
    z0 = hfun(x0, U[:,1], D[:,1], p)
    nz = length(z0)
    Z = zeros(nz, length(t_vec)-1)
    
    for i in 1:(length(t_vec)-1)
        f = ffun(X[:,i], U[:,i], D[:,i], p)
        X[:,i+1] = X[:,i] + dt_sim * f
        Z[:,i] = hfun(X[:,i], U[:,i], D[:,i], p)
    end
    return X, Z
end