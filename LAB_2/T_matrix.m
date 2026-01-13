function T = T_matrix(t)
    % t is a column vector of 16 real variables
    T = [
        t(1), t(5) + 1i*t(6), t(7) + 1i*t(8), t(9) + 1i*t(10);
        0, t(2), t(11) + 1i*t(12), t(13) + 1i*t(14);
        0, 0, t(3), t(15) + 1i*t(16);
        0, 0, 0, t(4)
    ];
end