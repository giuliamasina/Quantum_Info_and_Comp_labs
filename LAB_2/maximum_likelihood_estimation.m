function rho_mle_matrices = maximum_likelihood_estimation(N_values, N_names, gammas)

    num_sets = size(N_values, 2); % Number of columns is the number of data sets
    rho_mle_matrices = cell(num_sets, 1);
    disp('Starting Quantum Maximum Likelihood Estimation...');
    disp('-----------------------------------------------');

    for k = 1:num_sets

        N = N_values(:, k); 
        N_name = N_names{k};
        
        fprintf('Processing: **%s**\n', N_name);
        
        % optimization
        t_optimum = optimize_variables(N, gammas);
        
        % compute final density matrix
        rho_mle = compute_rho_mle(t_optimum);
        
        % check physicality
        if is_density_matrix(rho_mle)
            disp('The density matrix after MLE is **physical** (positive semi-definite and trace 1).');
        else
            disp('The density matrix after MLE is **NOT physical**. (Consider using fmincon for constraints).');
        end
        
        rho_mle_matrices{k} = rho_mle;
        disp('rho_mle =');
        disp(rho_mle);
        fprintf('\n');
    end

    

end