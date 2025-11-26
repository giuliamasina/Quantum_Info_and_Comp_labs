function vonneumann_entropy(rho_mle_matrices)

    disp('Computing Von-Neumann Entropy:');
    disp('------------------------------');

    for i = 1:length(rho_mle_matrices)
        
        vn_entropy_value = compute_entropy(rho_mle_matrices{i});
        fprintf('Von-Neumann entropy (rho_mle_%d): %f\n', i-1, vn_entropy_value);

    end

end

function entropy = compute_entropy(rho)

    %COMPUTE_VON_NEUMANN_ENTROPY Computes the Von-Neumann entropy of a density matrix rho.
    %   S(rho) = -sum(lambda_i * log2(lambda_i))
    %   where lambda_i are the eigenvalues of rho.

    % 1. Compute eigenvalues of the Hermitian matrix rho
    % 'eig' computes eigenvalues. For a Hermitian matrix (like a density matrix),
    % it returns real eigenvalues.
    eigenvalues = eig(rho);
    
    % 2. Find non-zero eigenvalues (important for the log, as log(0) is -Inf)
    % The tolerance should be small but greater than machine precision (eps).
    % Use a small tolerance instead of strictly '!= 0' to account for numerical errors.
    %tolerance = 1e-12;
    %non_zero_eigenvalues = eigenvalues(abs(eigenvalues) > tolerance);
    non_zero_eigenvalues = eigenvalues(abs(eigenvalues) ~= 0);
    
    % 3. Compute the term -lambda_i * log2(lambda_i) for each non-zero eigenvalue
    % 'log2' is the base-2 logarithm in MATLAB. The '.' operator performs element-wise multiplication.
    terms = -non_zero_eigenvalues .* log2(non_zero_eigenvalues);
    
    % 4. Sum the terms to get the entropy
    entropy = sum(terms);

end