%First to execute

% starting values
I0 = max(WaveplateA.coincidences) - min(WaveplateA.coincidences);
[~, idx] = min(WaveplateA.coincidences);
theta0 = WaveplateA.alpha(idx);
offset = min(WaveplateA.coincidences)

% Define custom fit type (Malus' Law)
ft = fittype('I0 * cosd(2*(x - theta0)).^2 + offset','independent', 'x');
startingValues = [I0, theta0, offset];
fittedModel = fit(WaveplateA.angles, WaveplateA.coincidences, ft, StartPoint=startingValues);

angle_range = linspace(min(WaveplateA.angles),max(WaveplateA.angles),2000);
modelEval = feval(fittedModel, angle_range);
[min_sin, index_min_sin] = min(modelEval);
angle_min = angle_range(index_min_sin);

figure;
plot(fittedModel)
hold on
plot(WaveplateA.angles,WaveplateA.coincidences,"bo")
plot(angle_min, min_sin, 'r.', 'MarkerSize', 20);
ylim([min(modelEval)-500, max(modelEval) + 500])
xlim([min(WaveplateA.angles), max(WaveplateA.angles)])
title('WPA calibration');
xlabel('angle [Deg]');
ylabel('coincidences');
legend('Data', 'Fitted function');
hold off

