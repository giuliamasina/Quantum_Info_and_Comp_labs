function concurrence(rho_mle_matrices)

    disp('Computing Concurrence:');
    disp('----------------------');

    for i = 1:length(rho_mle_matrices)
        
        concurrence_value = compute_concurrence(rho_mle_matrices{i});
        fprintf('Concurrence (rho_mle_%d): %f\n', i-1, concurrence_value);

    end

end