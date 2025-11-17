function rho = compute_rho_lin_inv(N, psi_states, gammas)
% COMPUTE_RHO_LIN_INV Determines the density matrix (rho) by linear inversion.
%
%   rho = COMPUTE_RHO_LIN_INV(N, psi_states, gammas)
%
%   Inputs:
%       N          - A 16-element column vector containing the measured counts 
%                    for the 16 configurations.
%       psi_states - A 4x16 matrix where each column is a state vector |psi_i>.
%       gammas     - A 1x16 cell array containing the 4x4 gamma matrices.
%
%   Output:
%       rho        - The 4x4 density matrix.

    % N should be a 16x1 column vector for matrix multiplication later
    if size(N, 1) < size(N, 2)
        N = N.';
    end

    % Check the dimensions
    if size(psi_states, 2) ~= 16 || length(gammas) ~= 16 || length(N) ~= 16
        error('Input dimensions must be 16 for psi_states columns, gammas cells, and N elements.');
    end
    
    % --- 1. Compute the B matrix ---
    
    % Initialize B (16x16 complex matrix)
    B = zeros(16, 16);
    
    for i = 1:16
        for j = 1:16
            % psi_states(:, i) is the column vector |psi_i>
            psi_i = psi_states(:, i);
            % gammas{j} is the operator Gamma_j
            gamma_j = gammas{j};
            
            % Python's np.conj(psi_states[i].T) is MATLAB's psi_i' (Hermitian conjugate)
            % Python's np.dot(A, B) is MATLAB's A * B
            
            % B[i, j] = <psi_i | Gamma_j | psi_i>
            % Note: MATLAB's index i is 1-based, Python's is 0-based, but 
            % the loops map directly 1:16.
            B(i, j) = psi_i' * gamma_j * psi_i;
        end
    end
    
    % --- 2. Compute the inverse of B (B_inv) ---
    B_inv = inv(B);
    
    % --- 3. Compute the M matrices (M_i) ---
    
    % Initialize M (a 4x4x16 3D array, where M(:, :, i) is the matrix M_i)
    M = zeros(4, 4, 16);
    
    % In MATLAB, it's generally cleaner to loop over the index 'mu' first.
    for mu = 1:16
        M_mu = zeros(4, 4); % Initialize the current M matrix
        for nu = 1:16
            % The Python code uses B_inv[nu, mu] (row nu, column mu)
            scalar_coeff = B_inv(nu, mu);
            gamma_nu = gammas{nu};
            
            % M[mu] += B_inv[nu, mu] * gammas[nu]
            M_mu = M_mu + scalar_coeff * gamma_nu;
        end
        M(:, :, mu) = M_mu; % Store the resulting matrix in the 3D array
    end
    
    % --- 4. Compute the density matrix rho ---
    
    % Initialize rho (4x4 complex matrix)
    rho = zeros(4, 4);
    
    % rho = sum(M[i] * N[i])
    for i = 1:16
        % M(:, :, i) is the matrix M_i
        % N(i) is the scalar count N_i
        rho = rho + M(:, :, i) * N(i);
    end
    
    % --- 5. Normalization ---
    
    % The Python code normalizes by the sum of the first 4 counts: 
    % rho /= np.sum(N[:4])
    % This is likely a typo in the original Python for a general QST normalization.
    % However, to be strictly faithful to the provided code:
    
    normalization_factor = sum(N(1:4));
    
    % If the input N represents measured counts, they should be real and non-negative.
    rho = rho / normalization_factor; 
    
    % A more standard normalization for the density matrix is trace(rho) = 1.
    % If you wanted the standard normalization:
    % rho = rho / trace(rho);

end