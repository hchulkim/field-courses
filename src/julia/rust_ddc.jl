using DataFrames, CSV, LinearAlgebra, Statistics, SparseArrays

## q1

# set parameter values
theta1 = 0.0049
theta2 = 0.0002
theta3 = 0.35
R = 2.202
beta = 0.95
gamma = 0.5772156649
K = 90
tol = 1e-12
maxit = 20000

true_theta = [theta1; theta2; R]

x = collect(0:K-1)
c = theta1 .* x .+ theta2 .* (x .^ 2)  

# transition: F0
F0 = spdiagm(
    0 => theta3 .* ones(K),
    1 => (1 - theta3) .* ones(K-1),
)
F0[end, end] = 1.0

# transition: F1
F1 = spzeros(K, K)
F1[:, 1] .= theta3
if K >= 2
    F1[:,2] .= 1 - theta3
end

# Elementwise, numerically-stable log sum exp
ev_logsum(a, b) = begin
    m =  log.(exp.(a) + exp.(b))
end

# q1
EV = zeros(K, 1)
it = 0
dist = Inf

while dist > tol && it < maxit

    it += 1
    delta0 = .-c .+ beta .* (F0 * EV)
    delta1 = .-R .+ beta .* (F1 * EV)
    EV_new = ev_logsum(delta0, delta1)
    dist = norm(EV_new - EV)
    EV = EV_new
end

detal0 = .-c .+ beta .* (F0 * EV)
detal1 = .-R .+ beta .* (F1 * EV)
EV_vfi = copy(EV)

if it == maxit
    @warn "Maximum number of iterations reached"
end

# save EV_vfi
CSV.write("input/temp/EV_vfi.csv", DataFrame(EV_vfi, :auto))

# Optimal policy
denom = exp.(detal0) .+ exp.(detal1)

P1_vfi = exp.(detal1) ./ denom
P0_vfi = exp.(detal0) ./ denom


## q2
it = 0

P0 = 0.5 .* ones(K, 1)
P1 = 0.5 .* ones(K, 1)

Id = spdiagm(
    0 => ones(K)
)

D0 = spdiagm(
    0 => 0.5 .* ones(K)
)

D1 = spdiagm(
    0 => 0.5 .* ones(K)
)



while dist > tol && it < maxit

    it += 1
    ev_new = Inverse(Id - beta * (D0 * F0 + D1 * F1)) * ()
end