% Cleaning: load from the results directory to postprocess a different file
close all; clear; clc

%% Error Propagation
addpath('LiteratureModels')

D = linspace(0.001, 0.02, 20);
n_samp = 100;

% Covariance Matrix
M = [phi1.var, phi1.cov; phi2.cov, phi2.var];

% Vector of means
means = [phi1.best; phi2.best];

% Cholesky decomposition
R = chol(M);

% Monte Carlo sampling
vt_MC = nan(length(D), n_samp);
Y_MC = nan(length(D), size(M, 1), n_samp);

wb = waitbar(0, 'Generating correlated random samples.. ');
for i = 1:1:length(D)
    % generate uncorrelated random sequence
    x = rand(size(M, 1), n_samp);

    % correlate the sequence
    y = means + R' * x;

    % generate the correlated sample
    for j = 1:1:n_samp
        vt_MC(i, j) = vt(D(i), y(1, j), y(2, j), model);
    end
    
    % storing correlated random variables
    Y_MC(i, :, :) = y;
    waitbar(i/length(D), wb)
end
close(wb)

% Statistics calculation
vt_mean = mean(vt_MC, 2);
vt_std = std(vt_MC, 0, 2);

% result
dv = NaN*ones(1, length(data));
vt_ref = NaN*ones(1, length(data));
for i = 1:1:length(data)
    vt_ref(i) = data(i).vt;
    dv(i) = data(i).dv;
end

figure()
plot(dv, vt_ref, 'k')
hold on
errorbar(D, vt_mean, vt_std);
legend(data(1).name, model)

figure()
surf(X, Y, P{1})
view([0 90])
hold on
scatter3(phi1.best(1), phi2.best(1), max(max(P{1})), 'r', 'filled')
scatter3(Y_MC(1, 1, :), Y_MC(1, 2, :), ones(1, n_samp).*max(max(P{1})), 'g')