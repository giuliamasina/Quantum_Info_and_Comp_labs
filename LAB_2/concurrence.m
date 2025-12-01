function concurrence(rho_mle_matrices)

    disp('Computing Concurrence:');
    disp('----------------------');

    for i = 1:length(rho_mle_matrices)
        
        concurrence_value = compute_concurrence(rho_mle_matrices{i});
        fprintf('Concurrence (rho_mle_%d): %f\n', i-1, concurrence_value);

    end

end

function concur = compute_concurrence(rho)

    sigma = [0, 0, 0, -1; 
             0, 0, 1, 0; 
             0, 1, 0, 0; 
             -1, 0, 0, 0];
    rho_tilde = sigma * conj(rho) * sigma;
    sqrtm_rho = sqrtm(rho);
    M = sqrtm_rho * rho_tilde * sqrtm_rho;
    omega_hat = sqrtm(M);
    eigenvalues_full = eig(omega_hat);
    eigenvalues = real(eigenvalues_full);
    eigenvalues(eigenvalues < 0) = 0; % Set negative to zero
    eigenvalues = eigenvalues / sum(eigenvalues); % Non-standard normalization
    eigenvalues = sort(eigenvalues, 'descend'); % Sort descending
    concur = max(0, eigenvalues(1) - eigenvalues(2) - eigenvalues(3) - eigenvalues(4));

end