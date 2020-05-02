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
save('HolzerSommerfeld/All_001-1_05-1_200s_x2.mat')

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
save(['HoltzerSommerfeld/Seq' num2str(step) '_20_001-1_05-1_50s.mat'])
% save('test/test.mat')