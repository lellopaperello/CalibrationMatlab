% Calibration of snowflakes shape parameters using Bayesian approach.
% Double parameter estimation. 
% Governing equation: equilibrium of a point mass in free fall
%       1/2 rho_a vt^2 S(dv, Phi) cD(Re(vt, dv), model, Phi) = 
%       V(dv) (rho - rho_a) g 
% Formulation: vt = f(Phi, model, dv)
% vt = terminal velocity (experimental data)
% model = Ganser (one of the better fitting models)
% dv = diameter of the volume-equivalent sphere (experimental data)
% Phi = Sphericity (A_eq sphere / A_particle) ~ free parameter

%% All data at once: searching for a single best estimation of (Phi, Ar)
%  One-step data analysis
% OBS: A single couple (Phi, Ar) isn't enough to properly describe vt in
%      the whole dv range => overdetermined problem
close all; clear; clc %#ok<*NOPTS>
input_data

P = cell(1, 1);
P{1} = posterior(phi1.vec, phi2.vec, model, data);
save('results/All_001-03_001-02_100s_x2.mat')

%% Single datum analysis: searching for the best estimation of (Phi, Ar)_i
%  as a function of dv_i
% OBS: A whole set of (Phi, Ar)_i is optimal for each single dv_i
%      => underdetermined problem
close all; clear; clc %#ok<*NOPTS>
input_data

P = cell(1, length(data));
for i = 1:1:length(data)
    P{i} = posterior(phi1.vec, phi2.vec, model, data(i));
end
save('results/Seq10_001-1_001-1_100s.mat')

%% Data by small groups: searching for the best estimation of (Phi, Ar)_i
%  as a piecewise function of dv_i
close all; clear; clc %#ok<*NOPTS>
input_data

step = 3;

N = floor((length(data)-1) / step);
P = cell(1, N);
wb = waitbar(0, 'Passing data samples.. ');
pos = get(wb, 'position');
pos = [pos(1) pos(2) + 1.5*pos(4) pos(3) pos(4)];
set(wb, 'position', pos, 'doublebuffer', 'on')
for i = 1:1:N
    index = ((i-1)*step + 1):1:i*step + 1;
    P{i} = posterior(phi1.vec, phi2.vec, model, data(index));
    waitbar(i/N, wb)
end
close(wb)
save(['results/Seq' num2str(step) '_20_001-1_001-1_100s.mat'])
% save('test/test.mat')

%% Cleaning: load from "results" directory to postprocess a different file
close all; clear; clc
%% Postprocessing

phi1.EV = zeros(1, length(P));
phi2.EV = zeros(1, length(P));
for i = 1:1:length(P)
    % maximum/a
    [imax, jmax] = find(P{i} == max(max(P{i})));
    phi1.EV(i) = phi1.vec(jmax);
    phi2.EV(i) = phi2.vec(imax);

    % Posterior
    figure()
    [X, Y] = meshgrid(phi1.vec, phi2.vec);
    surf(X, Y, P{i})
    hold on
    scatter3(phi1.EV(i), phi2.EV(i), max(max(P{i})), 'r', 'filled')

    title('Posterior')
    xlabel(phi1.name)
    ylabel(phi2.name)
%     view([0 90])
    c = colorbar;
    c.Label.String = ['$ P(' phi1.name ', ' phi2.name ') | \{data\}) $'];
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 12;
end

% Logarithm of the posterior
% figure()
% [X, Y] = meshgrid(phi1.vec, phi2.vec);
% surf(X, Y, logP)
% hold on
% scatter(phi2.EV, phi1.EV)
% 
% title('Posterior')
% xlabel(phi1.name)
% ylabel(phi2.name)
% view([0 90])
% axis([phi1.min phi1.max phi2.min phi2.max -5e3 max(max(logP))])
% caxis([-5e3 max(max(logP))]);
% c = colorbar;
% c.Label.String = ['$ ln(P(' phi1.name ', ' phi2.name ') | \{data\})) $'];
% c.Label.Interpreter = 'latex';
% c.Label.FontSize = 12;

%% Result
vt_best1 = cell(1, length(P));
dv = NaN*ones(1, length(data));
vt_ref = NaN*ones(1, length(data));
for i = 1:1:length(data)
    vt_ref(i) = data(i).vt;
    dv(i) = data(i).dv;
    for j = 1:1:length(P)
        vt_best1{j}(i) = vt(data(i).dv, phi1.EV(j), phi2.EV(j), model);
    end
end

figure()
Legend = cell(1, length(P) + 1);
plot(1e3*dv, vt_ref, '--k')
Legend{1} = data(1).name;
hold on
for j = 1:1:length(P)
    plot(1e3*dv, vt_best1{j})
    Legend{j + 1} = [model ', d_v = ' num2str(data(j).dv)];
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
    plot(dv, phi1.EV, dv, phi2.EV)
    xlabel('d_v [mm]')
    legend(phi1.name, phi2.name)
end
