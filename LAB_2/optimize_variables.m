function t_opt = optimize_variables(N, gammas)
    % Initial guess for the 16 variables (t1 to t16).
    initial_guess = randn(16, 1); 
    
    % The objective function for fminsearch is the likelihood (Chi-Squared), 
    % which we want to MINIMIZE.
    objective_fun = @(t) -likelihood(t, N, gammas);

    options = optimset('MaxFunEvals', 2e5, 'MaxIter', 2e5);

    % Perform the optimization using the Nelder-Mead method
    t_opt = fminsearch(objective_fun, initial_guess, options);
end