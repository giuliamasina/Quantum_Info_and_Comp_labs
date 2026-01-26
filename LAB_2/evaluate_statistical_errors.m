function evaluate_statistical_errors(rho_mle_matrices, N_values, N_names, gammas)

    num_simulations = 150;
    num_N = length(N_names);         
    %num_sets = length(rho_mle_matrices);  
    
    simulated_fidelities   = zeros(num_simulations, num_N);
    simulated_concurrences = zeros(num_simulations, num_N);
    
    for sim_idx = 1:num_simulations
            if mod(sim_idx, 10) == 0
                fprintf("Progress: %d", sim_idx);
                fprintf('\n');
            end
            
            for i = 1:num_N
                % mean
                N_mean = N_values(:, i);
                simulated_counts = poissrnd(N_mean);
                
                % MLE
                t_optimum   = optimize_variables(simulated_counts, gammas);
                rho_mle_sim = compute_rho_mle(t_optimum);
                
                %fidelities
                rho_target   = rho_mle_matrices{i};
                fidelity_mle = fidelity_for_statistical_errors(rho_mle_sim, rho_target);
                simulated_fidelities(sim_idx, i) = fidelity_mle;
                
                % concurrences
                concurrence_mle = compute_concurrence(rho_mle_sim);
                simulated_concurrences(sim_idx, i) = concurrence_mle;
            end
    end
    

    % calculate the statistical errors as the standard deviations of the simulated results"
    fidelity_errors    = std(simulated_fidelities, 0, 1);    
    concurrence_errors = std(simulated_concurrences, 0, 1);
    
    fprintf('=== Simulation Results ===\n');
    for N_index = 1:num_N
        current_name = N_names{N_index};
        
        fprintf('Dataset: %s\n', current_name);
        
        fprintf('  Fidelity errors:   ');
        fprintf('%.5f  ', fidelity_errors(N_index));
        fprintf('\n');
        
        fprintf('  Concurrence errors:');
        fprintf('%.5f  ', concurrence_errors(N_index));
        fprintf('\n\n');
    end



end