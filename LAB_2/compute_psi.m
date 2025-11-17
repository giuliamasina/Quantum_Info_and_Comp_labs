
function psi = compute_psi(h_1, q_1, h_2, q_2)
% compute_psi computes the |psi> state vector based on waveplate angles.

    % Complex unit is 1i in MATLAB
    a_1 = (1 / sqrt(2)) * (sin(2 * h_1) - 1i * sin(2 * (h_1 - q_1)));
    b_1 = (-1 / sqrt(2)) * (cos(2 * h_1) + 1i * cos(2 * (h_1 - q_1)));
    
    a_2 = (1 / sqrt(2)) * (sin(2 * h_2) - 1i * sin(2 * (h_2 - q_2)));
    b_2 = (-1 / sqrt(2)) * (cos(2 * h_2) + 1i * cos(2 * (h_2 - q_2)));
    
    % Compute the four coefficients
    first_term = a_1 * a_2;
    second_term = a_1 * b_2;
    third_term = b_1 * a_2;
    fourth_term = b_1 * b_2;
    
    % Construct the 4x1 state vector (column vector)
    % MATLAB prefers column vectors for state representations
    psi = [first_term; second_term; third_term; fourth_term];
    
    % If you need a row vector, use:
    % psi = [first_term, second_term, third_term, fourth_term];
end