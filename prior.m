function P = prior(phi1, phi2, type)
input.m
switch type
    case 'uniform'
        P = zeros(length(phi1), length(phi2));
        constant = 1 / ((phi1_M - phi1_m) * (phi2_M - phi2_m));
        P( (phi1 > phi1_m && phi1 < phi1_M) && ...
           (phi2 > phi2_m && phi2 < phi2_M) ) = constant;
end