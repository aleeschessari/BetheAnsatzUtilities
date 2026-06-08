module LiebLiniger

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

# TODO: Add Lieb-Liniger model specific implementations
# These functions should mirror the Heisenberg_XXX interface but implement
# the quantum gas dynamics instead of spin chain physics

function create_available_array(N, r)
    error("LiebLiniger.create_available_array not yet implemented")
end

function get_ground_state(N; tol=1e-12, maxiter=1000)
    error("LiebLiniger.get_ground_state not yet implemented")
end

function range_I(N, r)
    error("LiebLiniger.range_I not yet implemented")
end

function get_energy_triplet(I_array_available, r, N, x; tol=1e-12, maxiter=1000)
    error("LiebLiniger.get_energy_triplet not yet implemented")
end

function get_momentum_index_triplet(I_array_available, r, N, x)
    error("LiebLiniger.get_momentum_index_triplet not yet implemented")
end

function energy(z, N)
    error("LiebLiniger.energy not yet implemented")
end

function phi(x)
    error("LiebLiniger.phi not yet implemented")
end

function iterative_procedure(I_array, r, N; tol=1e-12, maxiter=10000)
    error("LiebLiniger.iterative_procedure not yet implemented")
end

function d(x, N)
    error("LiebLiniger.d not yet implemented")
end

function K_function(z)
    error("LiebLiniger.K_function not yet implemented")
end

function kappa_function(z)
    error("LiebLiniger.kappa_function not yet implemented")
end

function K_matrix(z, N)
    error("LiebLiniger.K_matrix not yet implemented")
end

function L_cursive(z)
    error("LiebLiniger.L_cursive not yet implemented")
end

function log_L_cursive(z)
    error("LiebLiniger.log_L_cursive not yet implemented")
end

function kappa_cursive(z)
    error("LiebLiniger.kappa_cursive not yet implemented")
end

function log_kappa_cursive(z)
    error("LiebLiniger.log_kappa_cursive not yet implemented")
end

function G_function(z)
    error("LiebLiniger.G_function not yet implemented")
end

function H_matrix(z0, z, N)
    error("LiebLiniger.H_matrix not yet implemented")
end

function H_matrix_simplified_minus(z0, z, N)
    error("LiebLiniger.H_matrix_simplified_minus not yet implemented")
end

function H_matrix_simplified_plus(z0, z, N)
    error("LiebLiniger.H_matrix_simplified_plus not yet implemented")
end

function form_factor_simplified_minus(z0, z, N)
    error("LiebLiniger.form_factor_simplified_minus not yet implemented")
end

end # module LiebLiniger
