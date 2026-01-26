function rho_mle = compute_rho_mle(t)
    T = T_matrix(t);
    
    normalization_factor = sqrt(trace(T' * T));   
    %check for division by zero
    if normalization_factor < 1e-10
        rho_mle = eye(4)/4; 
        return;
    end
    
    % normalize T 
    T_norm = T / normalization_factor;
    
    % compute rho
    rho_mle = T_norm' * T_norm; 
end