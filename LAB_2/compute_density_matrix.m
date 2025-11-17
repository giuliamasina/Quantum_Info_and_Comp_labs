
coinc_phi_p = Cartel2.coinc_phi_p;
coinc_phi_m = Cartel2.coinc_phi_m;
coinc_dec = Cartel2.coinc_dec;

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

% Define gamma matrices (tensor product of Pauli matrices)

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

% --- Collect the matrices into a cell array ---

% In MATLAB, a "list" of matrices/arrays (where all elements aren't the same size) 
gammas = {gamma_1, gamma_2, gamma_3, gamma_4, gamma_5, gamma_6, ...
          gamma_7, gamma_8, gamma_9, gamma_10, gamma_11, gamma_12, ...
          gamma_13, gamma_14, gamma_15, gamma_16};

% Accessing a matrix: use curly braces, e.g., gammas{1}
% Accessing an element: use gammas{1}(row, col)

% Get the number of configurations (rows)
num_configurations = size(waveplates_configuration, 1);

% Pre-allocate psi_states matrix (4 rows, num_configurations columns)
% Each column will hold one 4-element psi state vector.
psi_states = complex(zeros(4, num_configurations));

% Loop through each row of the configuration matrix
for nu = 1:num_configurations
    % Extract the four angle variables for the current configuration
    angles = waveplates_configuration(nu, :);
    h1 = angles(1);
    q1 = angles(2);
    h2 = angles(3);
    q2 = angles(4);
    
    % Compute the psi state using the function
    psi = compute_psi(h1, q1, h2, q2);
    
    % Store the resulting 4x1 vector as a column in psi_states
    psi_states(:, nu) = psi; 
end

num_sets = size(N_values, 3); % Number of columns is the number of data sets
rho_lin_inv_matrices = cell(1, num_sets);

for k = 1:num_sets
    % Extract the current N vector (column k)
    N = N_values(:, k); 
    
    % Get the current name
    N_name = N_names{k};
    
    % 1. Determine density matrix by linear inversion
    rho_lin_inv = compute_rho_lin_inv(N, psi_states, gammas);
    
    % 2. Check for physicality of the density matrix
    is_physical = is_density_matrix(rho_lin_inv);
    
    % 3. Display the result
    if is_physical
        fprintf('The density matrix for %s is physical.\n', N_name);
    else
        fprintf('The density matrix for %s is NOT physical.\n', N_name);
    end
    
    % 4. Store the matrix and print it
    rho_lin_inv_matrices{k} = rho_lin_inv;
    disp(rho_lin_inv);
    fprintf('\n'); % Print an empty line for separation
end

% The results are now stored in the cell array rho_lin_inv_matrices.

% To access the first calculated matrix:
% rho_A = rho_lin_inv_matrices{1};


