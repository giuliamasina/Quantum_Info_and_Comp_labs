function compute_fidelity(rho_mle_matrices)

    disp('Computing Fidelities:');
    disp('---------------------');

    for i = 1:length(rho_mle_matrices) 
            
            rho_1 = rho_mle_matrices{i};
            rho_2 = rho_mle_matrices{1};

            %compute the matrix square root of the first density matrix: sqrt(rho_1)
            sqrt_rho_1 = sqrtm(rho_1);
            
           % compute the argument inside the second square root: sqrt(rho_1) * rho_2 * sqrt(rho_1)
            M = sqrt_rho_1 * rho_2 * sqrt_rho_1;
            
            % compute the matrix square root of M
            sqrt_M = sqrtm(M);
            
            % compute the trace and take the real part 
            trace_value = real(trace(sqrt_M));
    
            % local fidelity function
            fidelity_value = trace_value^2;
            
            fprintf('Fidelity (rho_mle_%d, rho_mle_%d): %f\n', i, 1, fidelity_value);

    end


end