function coincidences = process_single_file(file)

    data = readtable(file);
    file_name = file(end-11:end);
    data.Properties.VariableNames = {'TimeTag', 'Channel'};
    %rescaling
    data.TimeTag = data.TimeTag * 81e-12;

    % i consider only the first second
    start_time = data.TimeTag(1);
    end_time = start_time + 1;
    data = data(data.TimeTag >= start_time & data.TimeTag <= end_time, :);

    %calculate time differences
    timetags = data.TimeTag;
    channels = data.Channel;

    delta_t = zeros(length(timetags),1);
    for i = 1:length(timetags) - 1
        delta_t(i) = (timetags(i+1)-timetags(i))*sign(channels(i+1)-channels(i));
    end
    delta_t = delta_t((delta_t ~= 0) & abs(delta_t) < 10e-9);

    % calculate histograms
    [counts, edges] = histcounts(delta_t, 80); 
    bin_centers = edges(1:end-1) + diff(edges)/2;
    x = bin_centers;
    y = counts;
    x = x(:);
    y = y(:);

    % preparing for gaussian fit (starting points and bounds)
    [max_counts, max_idx] = max(y);
    start_a1 = max_counts;          % amplitude
    start_b1 = x(max_idx);          % mean
    start_c1 = (max(x) - min(x)) / 10; % width
    lower_bounds = [0, min(x), 1e-12]; 
    upper_bounds = [max(y) * 2, max(x), max(x) - min(x)]; 
    
    % gaussian fit
    gaussian_fit = fit(x, y, 'gauss1', 'StartPoint', ...
        [start_a1, start_b1, start_c1],'Lower', lower_bounds, ...
        'Upper', upper_bounds);

    mean = gaussian_fit.b1;
    width = gaussian_fit.c1;
    sigma = width / sqrt(2);
    lower = mean - 3*sigma;
    upper = mean + 3*sigma;

    coincidences = sum(delta_t >= lower & delta_t <= upper);

    figure;
    bar(bin_centers, counts, 'hist');
    hold on;
    plot(gaussian_fit);
    xlabel('delta_t');
    ylabel('counts');
    title(file_name);
    hold off;
    % counts: number of delta_t values in each of the 80 bins, edges: length = 81, since there s one more edge than bins

end