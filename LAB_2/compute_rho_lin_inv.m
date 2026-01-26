function rho = compute_rho_lin_inv(N, psi_states, gammas)

    if size(N, 1) < size(N, 2)
        N = N.';
    end

    % check dimensions
    if size(psi_states, 2) ~= 16 || length(gammas) ~= 16 || length(N) ~= 16
        error('Input dimensions must be 16 for psi_states columns, gammas cells, and N elements.');
    end

    
    %normalization factor based on the sum of the first 4 counts
    norm_factor = sum(N(1:4));
   
    % normalization
    N = N / norm_factor;
    
    % initialize B (16x16 complex matrix)
    B = zeros(16, 16);
    
    for i = 1:16
        for j = 1:16        
            psi_i = psi_states(:, i);
            gamma_j = gammas{j};
            B(i, j) = psi_i' * gamma_j * psi_i;
        end
    end
    
    % inverse of B 
    B_inv = inv(B);
    
    % compute the M matrices (M_i)
    M = zeros(4, 4, 16);
    
    for mu = 1:16
        M_mu = zeros(4, 4);
        for nu = 1:16
            
            scalar_coeff = B_inv(nu, mu);
            gamma_nu = gammas{nu};
            
            M_mu = M_mu + scalar_coeff * gamma_nu;
        end
        M(:, :, mu) = M_mu; 
    end
    
    % compute the density matrix rho
    rho = zeros(4, 4);   
    for i = 1:16
        rho = rho + M(:, :, i) * N(i);
    end

end