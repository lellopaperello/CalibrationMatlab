% Input file for two-parameters data analysis
addpath('LiteratureModels')

% -------------------------------------------------------------------------
% AMBIENT PARAMETERS
% -------------------------------------------------------------------------
global rho_a mu g
rho_a = 1.225;      % air density               [ kg/m^3 ]
mu = 1.715e-5;      % air viscosity             [ kg/ms ]
g = 9.81;           % gravity field             [ N/kg ]

% -------------------------------------------------------------------------
% MODEL PARAMETERS
% -------------------------------------------------------------------------
model = 'Ganser';

% First parameter
phi1.name = '\Phi [-]';
phi1.min = 0.01;
phi1.max = 0.3;
phi1.vec = linspace(phi1.min, phi1.max, 200);

% Second parameter
phi2.name = 'd_n / d_v [-],';
phi2.min = 0.01;
phi2.max = 0.2;
phi2.vec = linspace(phi2.min, phi2.max, 200);

% -------------------------------------------------------------------------
% DATA
% -------------------------------------------------------------------------
Nsamples = 100;
% Structure declaration
structure.name = ' ';
structure.dv = 0;
structure.vt = 0;
structure.sigma = 0;

data = repmat(structure, 1, Nsamples);
dv = linspace(0.001, 0.02, Nsamples);
for i = 1:1:Nsamples
    data(i).name = 'Brandes - 2008, T = -10°C';
    data(i).dv = dv(i);
    data(i).vt = reference(data(i).dv, data(i).name);
    data(i).sigma = 2*0.17;
end

% for i = 1:1:Nsamples
%     data(i).name = 'Prova';
%     data(i).dv = dv(i);
%     data(i).vt = 0.6 * dv(i)^0.5;
%     data(i).sigma = 0.01;
% end

clearvars structure dv Nsamples