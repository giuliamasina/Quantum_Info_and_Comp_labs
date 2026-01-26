
coinc_phi_p = Lab2data.coinc_phi_p;
coinc_phi_m = Lab2data.coinc_phi_m;
coinc_dec = Lab2data.coinc_dec;

disp(Lab2data);

N_values = [coinc_phi_m, coinc_phi_p, coinc_dec];
N_names = {'coinc_phi_m','coinc_phi_p','coinc_dec'};

waveplates_configuration = [
    0, 0, 0, 0;
    0, 0, pi/4, 0;
    pi/4, 0, pi/4, 0;
    pi/4, 0, 0, 0;
    pi/4, pi/4, 0, 0;
    pi/4, pi/4, pi/4, 0;
    pi/8, pi/4, pi/4, 0;
    pi/8, pi/4, 0, 0;
    pi/8, pi/4, 0, -pi/4;
    pi/8, pi/4, pi/8, -pi/4;
    pi/4, pi/4, pi/8, -pi/4;
    0, 0, pi/8, -pi/4;
    pi/4, 0, pi/8, -pi/4;
    pi/4, 0, pi/4, -pi/4;
    0, 0, pi/4, -pi/4;
    pi/4, pi/4, pi/4, -pi/4
];

% gamma matrices (tensor product of Pauli matrices)

gamma_1 = 0.5 * [0, 1, 0, 0;
                 1, 0, 0, 0;
                 0, 0, 0, 1;
                 0, 0, 1, 0];

gamma_2 = 0.5 * [0, -1i, 0, 0;
                 1i, 0, 0, 0;
                 0, 0, 0, -1i;
                 0, 0, 1i, 0];

gamma_3 = 0.5 * [1, 0, 0, 0;
                 0, -1, 0, 0;
                 0, 0, 1, 0;
                 0, 0, 0, -1];

gamma_4 = 0.5 * [0, 0, 1, 0;
                 0, 0, 0, 1;
                 1, 0, 0, 0;
                 0, 1, 0, 0];

gamma_5 = 0.5 * [0, 0, 0, 1;
                 0, 0, 1, 0;
                 0, 1, 0, 0;
                 1, 0, 0, 0];

gamma_6 = 0.5 * [0, 0, 0, -1i;
                 0, 0, 1i, 0;
                 0, -1i, 0, 0;
                 1i, 0, 0, 0];

gamma_7 = 0.5 * [0, 0, 1, 0;
                 0, 0, 0, -1;
                 1, 0, 0, 0;
                 0, -1, 0, 0];

gamma_8 = 0.5 * [0, 0, -1i, 0;
                 0, 0, 0, -1i;
                 1i, 0, 0, 0;
                 0, 1i, 0, 0];

gamma_9 = 0.5 * [0, 0, 0, -1i;
                 0, 0, -1i, 0;
                 0, 1i, 0, 0;
                 1i, 0, 0, 0];

gamma_10 = 0.5 * [0, 0, 0, -1;
                  0, 0, 1, 0;
                  0, 1, 0, 0;
                  -1, 0, 0, 0];

gamma_11 = 0.5 * [0, 0, -1i, 0;
                  0, 0, 0, 1i;
                  1i, 0, 0, 0;
                  0, -1i, 0, 0];

gamma_12 = 0.5 * [1, 0, 0, 0;
                  0, 1, 0, 0;
                  0, 0, -1, 0;
                  0, 0, 0, -1];

gamma_13 = 0.5 * [0, 1, 0, 0;
                  1, 0, 0, 0;
                  0, 0, 0, -1;
                  0, 0, -1, 0];

gamma_14 = 0.5 * [0, -1i, 0, 0;
                  1i, 0, 0, 0;
                  0, 0, 0, 1i;
                  0, 0, -1i, 0];

gamma_15 = 0.5 * [1, 0, 0, 0;
                  0, -1, 0, 0;
                  0, 0, -1, 0;
                  0, 0, 0, 1];

gamma_16 = 0.5 * [1, 0, 0, 0;
                  0, 1, 0, 0;
                  0, 0, 1, 0;
                  0, 0, 0, 1];


gammas = {gamma_1, gamma_2, gamma_3, gamma_4, gamma_5, gamma_6, ...
          gamma_7, gamma_8, gamma_9, gamma_10, gamma_11, gamma_12, ...
          gamma_13, gamma_14, gamma_15, gamma_16};


% number of configurations (rows)
num_configurations = size(waveplates_configuration, 1);

psi_states = complex(zeros(4, num_configurations));

% loop through each row of the configuration matrix
for nu = 1:num_configurations
    % extract the four angle variables for the current configuration
    angles = waveplates_configuration(nu, :);
    h1 = angles(1);
    q1 = angles(2);
    h2 = angles(3);
    q2 = angles(4);
    
    % compute the psi state using the function
    psi = compute_psi(h1, q1, h2, q2);
    psi_states(:, nu) = psi; 
end

num_sets = size(N_values, 2);
rho_lin_inv_matrices = cell(1, num_sets);

for k = 1:num_sets
    
    N = N_values(:, k); 
    N_name = N_names{k};
    
    %determine density matrix by linear inversion
    rho_lin_inv = compute_rho_lin_inv(N, psi_states, gammas);
    
    % is it physical?
    is_physical = is_density_matrix(rho_lin_inv);
    
    if is_physical
        fprintf('The density matrix for %s is physical.\n', N_name);
    else
        fprintf('The density matrix for %s is NOT physical.\n', N_name);
    end
    
    rho_lin_inv_matrices{k} = rho_lin_inv;
    disp(rho_lin_inv);
    fprintf('\n'); 
end


% MAXIMUM-LIKELIHOOD ESTIMATION

rho_mle_matrices = maximum_likelihood_estimation(N_values, N_names, gammas);

% FIDELITY

compute_fidelity(rho_mle_matrices);

% VON-NEUMANN ENTROPY

vonneumann_entropy(rho_mle_matrices);

% CONCURRENCE

concurrence(rho_mle_matrices);

% STATISTICAL ERRORS

evaluate_statistical_errors(rho_mle_matrices, N_values, N_names, gammas);

