function fidelity_value = fidelity_for_statistical_errors(rho_1, rho_2)

            % compute the matrix square root of the first density matrix: sqrt(rho_1)
            sqrt_rho_1 = sqrtm(rho_1);
            
            % compute the argument inside the second square root: sqrt(rho_1) * rho_2 * sqrt(rho_1)
            M = sqrt_rho_1 * rho_2 * sqrt_rho_1;
            
            % compute the matrix square root of M
            sqrt_M = sqrtm(M);
            
            % compute the trace and take the real part (Trace(sqrt_M))
            trace_value = real(trace(sqrt_M));
    
            %local fidelity function
            fidelity_value = trace_value^2;

end