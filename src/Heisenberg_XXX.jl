module Heisenberg_XXX

export
    # Array/initialization functions
    create_available_array,
    get_ground_state,
    range_I,
    
    # Energy/momentum calculation
    get_energy_triplet,
    get_momentum_index_triplet,
    energy,
    
    # Iterative solver
    iterative_procedure,
    
    # Utility functions
    phi,
    d,
    
    # Kernel and matrix functions
    K_function,
    kappa_function,
    K_matrix,
    L_cursive,
    log_L_cursive,
    kappa_cursive,
    log_kappa_cursive,
    G_function,
    
    # H matrices
    H_matrix,
    H_matrix_simplified_minus,
    H_matrix_simplified_plus,
    
    # Form factors
    form_factor_simplified_minus

function create_available_array(N, r)
    I_max = N - r - 1
    I_min = -I_max

    start_val = I_min
    end_val = (I_max) 

    I_array_available = collect(I_min:2.0:I_max)/ 2.0

    return I_array_available
end

function get_ground_state(N; tol=1e-12, maxiter=1000)
    iseven(N) || throw(ArgumentError("N must be even, got N=$N"))
    r_ground = Int(N / 2)
    I_array_available_ground = create_available_array(N, r_ground)
    z_ground = iterative_procedure(I_array_available_ground, r_ground, N, tol=tol, maxiter=maxiter)
    E_ground = energy(z_ground, N)
    return E_ground
end

function range_I(N, r)
    two_I_max = Int(N - r - 1)
    return isodd(two_I_max) ? (two_I_max + 1) : (two_I_max + 1)
end

function get_energy_triplet(I_array_available, r, N, x; tol=1e-12, maxiter=1000)
    I_array_0 = I_array_available[I_array_available .!= I_array_available[x[1]]]
    I_array = I_array_0[I_array_0 .!= I_array_available[x[2]]]

    z = iterative_procedure(I_array, r, N, tol=tol, maxiter=maxiter)
    return energy(z, N)
end

function get_momentum_index_triplet(I_array_available, r, N, x)
    # momentum is 2 * pi / N * get_momentum_index
    I_array_0 = I_array_available[I_array_available .!= I_array_available[x[1]]]
    I_array = I_array_0[I_array_0 .!= I_array_available[x[2]]]

    return Int(mod(r * N / 2 - sum(I_array) - N^2 / 4, N))
end

function energy(z, N)
    return sum(-2 ./ (1 .+ z.^2)) / N
end

function phi(x)
    return 2 * atan(x)
end

function iterative_procedure(I_array, r, N; tol=1e-12, maxiter=10000)
    z = zeros(Float64, r)
    err = fill(Inf, r)

    iter = 0
    while maximum(abs.(err)) > tol && iter < maxiter
        iter += 1
        newz = similar(z)
        for i in 1:r
            iterable = [(z[i] - z[j]) / 2 for j in 1:r if j != i]
            # avoid calling a shadowed `phi` variable (could be a Matrix); inline the function
            newz[i] = tan(π / N * I_array[i] + 1 / (2 * N) * sum(2 .* atan.(iterable)))
            err[i] = (newz[i] - z[i])
        end
        z .= newz
    end

    if iter == maxiter
        @warn "iterative_procedure reached maxiter without converging" tol=tol maxiter=maxiter
    end

    return z
end

# d function
function d(x, N)
    return ((x - 1im) / (x + 1im))^N
end

# K_function
function K_function(z)
    return 4 / (4 + z^2)
end

# kappa_function
function kappa_function(z)
    return 2 / (1 + z^2)
end

# K_matrix
function K_matrix(z, N)
    K = complex(zeros(length(z), length(z)))
    for a in 1:length(z)
        for b in 1:length(z)
            if a != b
                K[a, b] = K_function(z[a] - z[b])
            else
                iterable = (K_function(z[a] - z[j]) for j in 1:length(z) if j != a)
                K[a, b] = N * kappa_function(z[a]) - sum(collect(iterable))
            end
        end
    end
    return K
end

# L_cursive
function L_cursive(z)
    prod = 1
    for i in 1:length(z)
        prod *= kappa_function(z[i])
    end
    return prod
end

# log L_cursive (stable for large systems)
function log_L_cursive(z)
    total = 0.0 + 0.0im
    for i in 1:length(z)
        total += log(kappa_function(z[i]))
    end
    return total
end

# kappa_cursive
function kappa_cursive(z)
    prod = 1
    for i in 1:length(z)
        for j in 1:length(z)
            if i < j
                prod *= K_function(z[i] - z[j])
            end
        end
    end
    return prod
end

# log kappa_cursive (stable for large systems)
function log_kappa_cursive(z)
    total = 0.0 + 0.0im
    for i in 1:length(z)
        for j in 1:length(z)
            if i < j
                total += log(K_function(z[i] - z[j]))
            end
        end
    end
    return total
end

# G_function
function G_function(z)
    return z / 2 + 1im
end

function H_matrix(z0, z, N)
    H = complex(zeros(length(z0), length(z0)))
    for a in 1:length(z0)
        for b in 1:length(z)
            prod0 = 1
            prod1 = 1
            for j in 1:length(z0)
                if j != a
                    prod0 *= (z0[j] - z[b] + 2*im)
                    prod1 *= (z0[j] - z[b] - 2*im)
                end
            end
            H[a, b] = 1im / (z0[a] - z[b]) * (prod0 - d(z[b], N) * prod1)
        end
    end
    return H
end

# H_matrix_simplified_minus
function H_matrix_simplified_minus(z0, z, N)
    H = complex(zeros(length(z0), length(z0)))
    for a in 1:length(z0)
        for b in 1:length(z)
            prod0 = 1
            prod1 = 1
            for j in 1:length(z0)
                if j != a
                    prod0 *= G_function(z0[j] - z[b])
                    prod1 *= conj(G_function(z0[j] - z[b]))
                end
            end
            H[a, b] = 1im / (z0[a] - z[b]) * (prod0 - d(z[b], N) * prod1)
        end
        H[a, end] = 1im * kappa_function(z0[a])
    end
    return H
end

# H_matrix_simplified_plus
function H_matrix_simplified_plus(z0, z, N)
    H = complex(zeros(length(z0), length(z0)))
    for a in 1:length(z0)
        for b in 1:length(z)
            prod0 = 1
            prod1 = 1
            for j in 1:length(z0)
                if j != a
                    prod0 *= G_function(z0[j] - z[b])
                    prod1 *= conj(G_function(z0[j] - z[b]))
                end
            end
            H[a, b] = 1im / (z0[a] - z[b]) * (prod0 - d(z[b], N) * prod1)
        end
        H[a, end] = 1im * kappa_function(z0[a])
    end
    return H
end

# form_factor_simplified_minus
function form_factor_simplified_minus(z0, z, N)
    #return N * (L_cursive(z) / L_cursive(z0)) * kappa_cursive(z0) * kappa_cursive(z) *
    #       abs(det(H_matrix_simplified_minus(z0, z, N)))^2 / (det(K_matrix(z0, N)) * det(K_matrix(z, N)))
    H = H_matrix_simplified_minus(z0, z, N)
    K0 = K_matrix(z0, N)
    K = K_matrix(z, N)

    log_prefactor = log(N) + log_L_cursive(z) - log_L_cursive(z0) +
                    log_kappa_cursive(z0) + log_kappa_cursive(z)

    return exp(log_prefactor + logdet(H) + logdet(conj.(H)) - logdet(K0) - logdet(K))
end

end # module Heisenberg_XXX
