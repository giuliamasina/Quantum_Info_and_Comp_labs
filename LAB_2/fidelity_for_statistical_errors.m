function fidelity_value = fidelity_for_statistical_errors(rho_1, rho_2)

            % 1. Compute the matrix square root of the first density matrix: sqrt(rho_1)
            sqrt_rho_1 = sqrtm(rho_1);
            
            % 2. Compute the argument inside the second square root: sqrt(rho_1) * rho_2 * sqrt(rho_1)
            % The matrix multiplication is done using the standard '*' operator
            M = sqrt_rho_1 * rho_2 * sqrt_rho_1;
            
            % 3. Compute the matrix square root of M
            sqrt_M = sqrtm(M);
            
            % 4. Compute the trace and take the real part (Trace(sqrt_M))
            % We use real() because the trace of a positive semi-definite matrix should be real,
            % but numerical noise can introduce a tiny imaginary part.
            trace_value = real(trace(sqrt_M));
    
            % Call the local fidelity function
            fidelity_value = trace_value^2;

end