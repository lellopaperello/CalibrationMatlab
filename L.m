function logP = L(phi1, phi2, model, data)
% Natural logarithm of the posterior. Since the likelihood function is a
% product of gaussians, its logarithm will be just a summation of quadratic
% terms. Since the logarithm is a monotonic function, the position of the
% maximum(a) will be inaltered.

logP = zeros(length(phi2), length(phi1));
for i = 1:1:length(phi2)
    for j = 1:1:length(phi1)
        for k = 1:1:length(data.dv)
            
            logP(i, j) = logP(i, j) - 0.5 * ((data.vt(k) ...
                       - vt(data.dv(k), phi1(j), phi2(i), model)) ...
                       / data.sigma(k))^2;
        end
    end
end