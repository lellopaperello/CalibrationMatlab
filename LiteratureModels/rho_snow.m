function rho = rho_snow(D)
% Brandes model for the calculation of the snowflakes density as function 
% of the ground temperature ~ Brandes et al. 2007, 
% "A Statistical and Physical Description of Hydrometeor Distributions
%  in Colorado Snowstorms Using a Video Disdrometer"

% D: equivalent volume diameter [ m ]
% rho: snowflake density [kg / m^3]

D = D*1e3;
rho = 0.178 * D.^(-0.922);
rho = rho*1e3;

end