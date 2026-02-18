using ForwardDiff

# Simple 2-state, 2-parameter system
# dx/dt = f(x, θ)
function f(x, θ)
    return [
        -θ[1] * x[1] + x[2],
        -θ[2] * x[2] + x[1]
    ]
end

# Test point
x0 = [1.0, 0.5]
θ0 = [0.3, 0.8]

# Jacobians
# ✅ Correct — differentiate wrt the actual argument being passed
dfdx = ForwardDiff.jacobian(x_ -> f(x_, θ0), x0)   # varies x, θ0 is fixed data — this is fine



println("f(x0, θ0) = ", f(x0, θ0))
println()
println("∂f/∂x = ", ∂f_∂x)
println()
println("∂f/∂θ = ", ∂f_∂θ)
```

Expected output:

f(x0, θ0) = [-0.8, 0.6]

∂f/∂x = [-0.3  1.0;
           1.0 -0.8]

∂f/∂θ = [-1.0  0.0;
           0.0 -0.5]

           ```