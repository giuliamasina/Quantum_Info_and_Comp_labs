function is_valid = is_density_matrix(rho)
% IS_DENSITY_MATRIX Checks if a given matrix (rho) satisfies the properties 
% of a physical density matrix.
%
%   Inputs:
%       rho - A square complex matrix.
%
%   Output:
%       is_valid - A logical value (true or false).

    % --- 1. Check Hermiticity (rho = rho') ---
    
    % rho' calculates the Hermitian conjugate (conjugate transpose) in MATLAB.
    rho_conj_transpose = rho';
    
    % allclose is replaced by abs(A - B) < tolerance. 
    % A common tolerance is set via eps (epsilon).
    tolerance = 1e-9; % Adjust tolerance as needed (1e-9 is usually safe)
    is_hermitian = all(all(abs(rho - rho_conj_transpose) < tolerance));
    
    % --- 2. Check Positivity (Eigenvalues >= 0) ---
    
    % eigvalsh is replaced by eig, specifically for Hermitian/Symmetric matrices.
    % We use 'hermitian' flag to ensure fast and stable calculation for Hermitian matrices.
    eigenvalues = eig(rho);
    
    % np.all(eigenvalues >= 0) is replaced by all(eigenvalues >= 0)
    %is_positive = all(eigenvalues >= -tolerance); % Check against a small negative tolerance
    is_positive = all(real(eigenvalues) >= -tolerance);
    
    % --- 3. Check Trace (Tr[rho] = 1) ---
    
    % trace(rho) computes the trace (sum of diagonal elements).
    rho_trace = trace(rho);
    
    % np.allclose(trace, 1) is replaced by abs(trace - 1) < tolerance
    is_trace_one = abs(rho_trace - 1) < tolerance;
    
    % --- 4. Final Check ---
    is_valid = is_hermitian && is_positive && is_trace_one;
    
    % You can also optionally display results for debugging:
    % fprintf('Hermitian Check: %d\n', is_hermitian);
    % fprintf('Positivity Check: %d\n', is_positive);
    % fprintf('Trace Check: %d\n', is_trace_one);

end