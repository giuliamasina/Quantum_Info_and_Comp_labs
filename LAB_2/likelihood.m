function chi_sq = likelihood(t, N, gammas)
    % N: Vector of measurement counts (experimental)
    % gammas: Cell array of measurement operators (M_m)
    
    rho_mle = compute_rho_mle(t);
    num_measurements = length(N);
    
    P = zeros(num_measurements, 1); % Predicted probabilities P_m = Tr(M_m * rho_mle)
    
    for m = 1:num_measurements
        % The predicted probability P_m is calculated as the trace
        P(m) = real(trace(gammas{m} * rho_mle)); % Use real() to handle numerical errors
    end
    
    % Compute the Chi-Squared value (Cost Function)
    % chi_sq = sum( (P - N)^2 / N )
    % Add a small epsilon (1e-10) to N to prevent division by zero, if needed.
    chi_sq = sum((P - N).^2 ./ (N + 1e-10)); 
end