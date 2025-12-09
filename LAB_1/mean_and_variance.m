function [mean, variance] = mean_and_variance(data, X, Y)

    reduced_data = data(data.X == X & data.Y == Y, :);
    N_total = sum(reduced_data.Coincidences);
    Npp = reduced_data(reduced_data.A == 0 & reduced_data.B == 0, :).Coincidences;
    Npm = reduced_data(reduced_data.A == 0 & reduced_data.B == 1, :).Coincidences;
    Nmp = reduced_data(reduced_data.A == 1 & reduced_data.B == 0, :).Coincidences;
    Nmm = reduced_data(reduced_data.A == 1 & reduced_data.B == 1, :).Coincidences;

    mean = (Npp + Nmm - Npm - Nmp)/N_total;
    term1_numerator = (N_total - mean)^2 * (Npp + Nmm);
    term2_numerator = (N_total + mean)^2 * (Npm + Nmp);
    denominator = N_total^4;
    
    % 3. Compute the variance
    variance = (term1_numerator + term2_numerator) / denominator;

end