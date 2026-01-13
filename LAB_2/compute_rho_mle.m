function rho_mle = compute_rho_mle(t)
    T = T_matrix(t);
    % T' in MATLAB is the Hermitian conjugate (T_dagger)
    T_dagger_T = T' * T; 
    % Normalization: rho_mle = T_dagger_T / trace(T_dagger_T)
    rho_mle = T_dagger_T / trace(T_dagger_T); 
end