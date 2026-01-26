function chi_sq = likelihood(t, N, gammas)
    % N: Vector of measurement counts (experimental)
    % gammas: Cell array of measurement operators (M_m)
    
    rho_mle = compute_rho_mle(t);
    num_measurements = length(N);

    % Calculate Total Counts (Normalization factor N from the text)
    total_counts = sum(N(1:4)); % Or sum(N) depending on your normalization definition
    
    P = zeros(num_measurements, 1); % Predicted probabilities P_m = Tr(M_m * rho_mle)
    
    for m = 1:num_measurements
        % The predicted probability P_m is calculated as the trace
        P(m) = real(trace(gammas{m} * rho_mle)); % Use real() to handle numerical errors
        %P(m) = total_counts * real(trace(psi_states(:, m) .* rho_mle)); % Use real() to handle numerical errors
        
        % Theoretical Probability
        %prob = real(trace(gammas{m} * rho_mle));
        
        % FIX 2: Convert Probability to Expected Counts (L_nu)
        % The text defines L_nu = N * <psi|rho|psi>
        %P(m) = total_counts * prob;
    end
    
    % Compute the Chi-Squared value (Cost Function)
    %chi_sq = sum( (P - N)^2 / N );
    % Add a small epsilon (1e-10) to N to prevent division by zero, if needed.
    chi_sq = sum((P - N).^2 ./ (N + 1e-10)); 
end