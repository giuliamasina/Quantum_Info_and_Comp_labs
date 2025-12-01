function rho_mle_matrices = maximum_likelihood_estimation(N_values, N_names, gammas)

    num_sets = size(N_values, 2); % Number of columns is the number of data sets
    rho_mle_matrices = cell(num_sets, 1);
    disp('Starting Quantum Maximum Likelihood Estimation...');
    disp('-----------------------------------------------');

    for k = 1:num_sets
        %N = N_values{k};
        %N_name = N_names{k};

        N = N_values(:, k); 
        N_name = N_names{k};
        
        fprintf('Processing: **%s**\n', N_name);
        
        % Perform Optimization
        t_optimum = optimize_variables(N, gammas);
        
        % Compute Final Density Matrix
        rho_mle = compute_rho_mle(t_optimum);
        
        % Check Physicality
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

% define the matrix T(t)
function T = T_matrix(t)
    % t is a column vector of 16 real variables
    T = [
        t(1), t(5) + 1i*t(6), t(7) + 1i*t(8), t(9) + 1i*t(10);
        0, t(2), t(11) + 1i*t(12), t(13) + 1i*t(14);
        0, 0, t(3), t(15) + 1i*t(16);
        0, 0, 0, t(4)
    ];
end

% define the function to compute rho_mle
function rho_mle = compute_rho_mle(t)
    T = T_matrix(t);
    % T' in MATLAB is the Hermitian conjugate (T_dagger)
    T_dagger_T = T' * T; 
    % Normalization: rho_mle = T_dagger_T / trace(T_dagger_T)
    rho_mle = T_dagger_T / trace(T_dagger_T); 
end

% define the function to compute the likelihood (Chi-Squared Cost Function)
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

% define the function to find the optimum set of variables (MLE)
function t_opt = optimize_variables(N, gammas)
    % Initial guess for the 16 variables (t1 to t16).
    initial_guess = randn(16, 1); 
    
    % The objective function for fminsearch is the likelihood (Chi-Squared), 
    % which we want to MINIMIZE.
    objective_fun = @(t) -likelihood(t, N, gammas);

    options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);

    % Perform the optimization using the Nelder-Mead method
    t_opt = fminsearch(objective_fun, initial_guess, options);
end