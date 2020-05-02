% Parametric studies of Holtzer and Sommerfeld cD models
close all; clear; clc
addpath('LiteratureModels')

% -------------------------------------------------------------------------
% AMBIENT PARAMETERS
% -------------------------------------------------------------------------
global rho_a mu g
rho_a = 1.225;      % air density               [ kg/m^3 ]
mu = 1.715e-5;      % air viscosity             [ kg/ms ]
g = 9.81;           % gravity field             [ N/kg ]

% Holtzer and Sommerfeld model: Parametric study over Phi and Phi_perp. 
% 
% vt = tree dimensional matrix storing the solution
% 1st dim: dv ~ indipendent variable
% 2nd dim: Phi ~ first free parameter
% 3rd dim: Phi_perp ~ second free parameter

% Indipendent variable
N1 = 1e2;
dv = linspace(0.001, 0.02, N1);

% Free parameters
N2 = 101;
N3 = 101;

% Physical boundaries
Phi_min = 0.01;
Phi_max = 1;
Phi_perp_min = 0.01;
Phi_perp_max = 1;

Phi = linspace(Phi_min, Phi_max, N2);
Phi_perp = linspace(Phi_perp_min, Phi_perp_max, N3);

% Solution computation
VT = zeros(N1, N2, N3);
wb = waitbar(0, 'Solution computation.. ');
for i = 1:1:N1
    for j = 1:1:N2
        for k = 1:1:N3
            VT(i, j, k) = vt(dv(i), Phi(j), Phi_perp(k), 'Holtzer&Sommerfeld');
        end
    end
    waitbar(i/N1)
end
close(wb)
save('results/VT_HoltzerSommerfeldStudy.mat')

%% Variable Phi @fixed Phi_perp
Phi_perp_values = [0.01 0.5 1];
step = 3;

kk = zeros(1, length(Phi_perp_values));
for i = 1:1:length(Phi_perp_values)
    [~, closest] = min(abs(Phi_perp -Phi_perp_values(i)));
    kk(i) = closest;
end


custom_color = jet(N2);
for k = kk
    figure()
    hold on
    for j = 1:step:N2
        plot(1e3*dv, squeeze(VT(:, j, k)), 'Color', custom_color(j, :))
    end
    colormap(jet(N2))
    cb = colorbar;
    title(['Phi_perp = ', num2str(Phi_perp(k))])
    xlabel('d_v [mm]')
    ylabel('v_t [m/s]')
    ylabel(cb, '\Phi [-]')
end

%% Variable Phi_perp @fixed Phi
Phi_values = [0.01 0.5 1];
step = 3;

jj = zeros(1, length(Phi_values));
for i = 1:1:length(Phi_values)
    [~, closest] = min(abs(Phi - Phi_values(i)));
    jj(i) = closest;
end

custom_color = jet(N3);
for j = jj
    figure()
    hold on
    for k = 1:1:N2
        plot(1e3*dv, squeeze(VT(:, j, k)), 'Color', custom_color(k, :))
    end
    colormap(jet(N2))
    cb = colorbar;
    title(['\Phi = ', num2str(Phi(j))])
    xlabel('d_v [mm]')
    ylabel('v_t [m/s]')
    ylabel(cb, '\Phi [-]')
end
