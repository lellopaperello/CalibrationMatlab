% Calibration of snowflakes shape parameters using Bayesian approach.
% Simple problem: single parameter estimation. 
% Governing equation: equilibrium of a point mass in free fall
%       1/2 rho_a vt^2 S(dv, Phi) cD(Re(vt, dv), model, Phi) = 
%       V(dv) (rho - rho_a) g 
% Formulation: vt = f(Phi, model, dv)
% vt = terminal velocity (experimental data)
% model = Chien (simple, single parameter model)
% cD = 30 / Re + 67.289*exp(-5.05*Phi)  ~  given
% dv = diameter of the volume-equivalent sphere (experimental data)
% Phi = Sphericity (A_eq sphere / A_particle) ~ free parameter

close all; clear; clc %#ok<*NOPTS>
addpath('LiteratureModels')

% -------------------------------------------------------------------------
% AMBIENT PARAMETERS
% -------------------------------------------------------------------------
rho_a = 1.225;      % air density               [ kg/m^3 ]
mu = 1.715e-5;      % air viscosity             [ kg/ms ]
g = 9.81;           % gravity field             [ N/kg ]

% Equation for terminal velocity (for Chien model is just quadratic)
A = 67.289;
B = 15*mu / rho_a;
C = g / (3*rho_a);
vt =@(Phi, dv) (-B./dv + sqrt(B^2./dv.^2 - A*exp(-5.05*Phi) *C.* Phi.* dv ...
                .* (rho_a - rho_snow(dv)))) ./ (A*exp(-5.05*Phi));

% References 
global Phi_min
global Phi_max
Phi_min = 0;
Phi_max = 2;

ref = 'Brandes - 2008, T = -10°C';
data.dv = linspace(0.001, 0.02, 1e2);
data.vt = reference(data.dv, ref);
data.sigma = 0.17;

% Likelihood ~ P( ({vt}, {dv}, {sigma}) | Phi) = P({data} | Phi) ----------
likelihood =@(Phi) prod(exp(-(data.vt - vt(Phi, data.dv)).^2 ./ 2*data.sigma.^2));

% Posterior ~ P(Phi | {data}) ---------------------------------------------
% Bayes' theorem: P(Phi | {data}) = P({data} | Phi) * P(Phi)
post =@(Phi) likelihood (Phi) * prior(Phi);

% One-step data analysis
phi = linspace(Phi_min, Phi_max, 1e3);
prob = ones(1, length(phi));
for i = 1:1:length(phi)
    prob(i) = post(phi(i));
end
% normalization condition ~ ind(prob(H)dH) = 1
prob = prob / trapz(phi, prob);

% Statistics
PHI = phi(prob == max(prob));

figure()
plot(phi, prob)
xline(PHI);
title('Final distribution')
xlabel('\Phi')
ylabel('pdf')
legend('posterior', '<\Phi>')

% Result
figure()
plot(1e3*data.dv, vt(PHI, data.dv))
hold on
plot(1e3*data.dv, data.vt)
title('Comparison')
xlabel('d_v [mm]')
ylabel('v_t [m/s]')
legend('v_t(\Phi_{best})', ref)

% Prior ~ P(Phi) ----------------------------------------------------------
function P = prior(Phi)
    global Phi_min
    global Phi_max
    if Phi >= Phi_min && Phi <= Phi_max
        P = 1 / (Phi_max - Phi_min);
    else
        P = 0;
    end
end