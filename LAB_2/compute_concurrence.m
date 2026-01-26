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
    eigenvalues(eigenvalues < 0) = 0; % set negative to zero
    eigenvalues = eigenvalues / sum(eigenvalues); % Non-standard normalization
    eigenvalues = sort(eigenvalues, 'descend'); % sort descending
    concur = max(0, eigenvalues(1) - eigenvalues(2) - eigenvalues(3) - eigenvalues(4));

end