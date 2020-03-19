function D = D_snow(T)
% Brandes model for the calculation of the equivalent volume diameter as
% function of the ground temperature ~ Brandes et al. 2007, 
% "A Statistical and Physical Description of Hydrometeor Distributions
%  in Colorado Snowstorms Using a Video Disdrometer"

% D: equivalent volume diameter [ m ]
% T: ground temperature [ K ]

D0 = -0.3475;
D1 = 0.13;
exp = 0.1733;
T0 = -49.04;

D = D1*(T-T0).^exp + D0;