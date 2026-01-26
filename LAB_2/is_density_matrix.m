function is_valid = is_density_matrix(rho)
    % check hermiticity (rho = rho')
    
    % conjugate transpose
    rho_conj_transpose = rho';
    tolerance = 1e-9; 
    is_hermitian = all(all(abs(rho - rho_conj_transpose) < tolerance));
    
    % check positivity (Eigenvalues >= 0) ---
    eigenvalues = eig(rho);
    is_positive = all(real(eigenvalues) >= -tolerance);
    
    % check trace equal to 1
    rho_trace = trace(rho);
    is_trace_one = abs(rho_trace - 1) < tolerance;
   
    
    is_valid = is_hermitian && is_positive && is_trace_one;
    

end