% Cleaning: load from the results directory to postprocess a different file
close all; clear; clc

%% Postprocessing
addpath('LiteratureModels')

% Statistcal parameters declaration
phi1.best = zeros(1, length(P));
phi2.best = zeros(1, length(P));

phi1.mean = zeros(1, length(P));
phi2.mean = zeros(1, length(P));

phi1.var = zeros(1, length(P));
phi2.var = zeros(1, length(P));

phi1.cov = zeros(1, length(P));
phi2.cov = zeros(1, length(P));

for i = 1:1:length(P)
    % integral
    I = trapz(phi2.vec, trapz(phi1.vec, P{i}, 2));
    
    % normalization
    if I ~= 1
        P{i} = P{i} / I;
    end
    
    % grid
    [X, Y] = meshgrid(phi1.vec, phi2.vec);

    % maximum
    [imax, jmax] = find(P{i} == max(max(P{i})));
    phi1.best(i) = phi1.vec(jmax);
    phi2.best(i) = phi2.vec(imax);

    % mean
    phi1.mean(i) = trapz(phi2.vec, trapz(phi1.vec, X.*P{i}, 2));
    phi2.mean(i) = trapz(phi1.vec, trapz(phi2.vec, Y.*P{i}), 2);
    
    % variance
    phi1.var(i) = trapz(phi2.vec, trapz(phi1.vec, ...
                        P{i}.*((X - phi1.mean(i)).^2), 2));
    phi2.var(i) = trapz(phi1.vec, trapz(phi2.vec, ...
                        P{i}.*((Y - phi2.mean(i)).^2), 2));
    
    % covariance
    phi1.cov(i) = trapz(phi2.vec, trapz(phi1.vec, ...
                        P{i}.*(X - phi1.mean(i)).*(Y - phi2.mean(i)), 2));
    phi2.cov(i) = phi1.var(i);
    
    % Posterior
    figure()
    surf(X, Y, P{i})
    hold on
    scatter3(phi1.best(i), phi2.best(i), max(max(P{i})), 'r', 'filled')
    xline(phi1.mean(i), '--r');
    yline(phi2.mean(i), '--r');
    
    title('Posterior')
    xlabel(phi1.name)
    ylabel(phi2.name)
    view([0 90])
    c = colorbar;
    c.Label.String = ['$ P(' phi1.name ', ' phi2.name ') | \{data\}) $'];
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 12;
end

% Display the results
disp('The statistics for phi1 are:')
disp(['Maximum = ' num2str(phi1.best)])
disp(['Mean = ' num2str(phi1.mean)])
disp(['Variance = ' num2str(phi1.var)])
disp(['Covariance = ' num2str(phi1.cov)])
disp(' ')
disp('The statistics for phi2 are:')
disp(['Maximum = ' num2str(phi2.best)])
disp(['Mean = ' num2str(phi2.mean)])
disp(['Variance = ' num2str(phi2.var)])
disp(['Covariance = ' num2str(phi2.cov)])

%% Result
vt_best1 = cell(1, length(P));
dv = NaN*ones(1, length(data));
vt_ref = NaN*ones(1, length(data));
for i = 1:1:length(data)
    vt_ref(i) = data(i).vt;
    dv(i) = data(i).dv;
    for j = 1:1:length(P)
        vt_best1{j}(i) = vt(data(i).dv, phi1.best(j), phi2.best(j), model);
    end
end

figure()
Legend = cell(1, length(P) + 1);
plot(1e3*dv, vt_ref, '--k')
Legend{1} = data(1).name;
hold on
for j = 1:1:length(P)
    plot(1e3*dv, vt_best1{j})
    Legend{j + 1} = model;
end
if length(P) > 1
    for j = 1:1:length(P)
        scatter(1e3*dv(j), vt_best1{j}(j))
    end
end

title('Comparison')
xlabel('d_v [mm]')
ylabel('v_t [m/s]')
legend(Legend)

if length(P) > 1
    figure()
    plot(phi1.best)
    hold on
    plot(phi2.best)
    xlabel('history')
    legend(phi1.name, phi2.name)
end
