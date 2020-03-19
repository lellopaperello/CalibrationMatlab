% "Is this a fair coin?"
% Implementation of the Example 1 (Ch. 2, pag. 14) of the book Data
% Analysis: A Bayesian Tutorial, D. S. Sivia and J. Skilling (2006)
close all; clear; clc %#ok<*NOPTS>

% Data samples [0 = tail  |  1 = head]
N = 1000; % number of tosses
w = 0.25; % bias weight: 0 ~ always TAIL, 1/2 ~ FAIR, 1 ~ always head

% Storing all the tosses
data = randsample([0 1], N, true, [1-w w]);

figure()
hist(data, [0 0.5 1])
axis([-1 2 0 N])
title('Data samples')
xlabel('value')
ylabel('number of events')

% Data Analysis

% Sequential approach: analyze the data as they arrive (each toss as 
% standalone). This is to avoid numerical errors (overshoots9 while 
% analyzing big numbers as a whole. Normalization at each analysis is
% required.

H = 0:0.001:1;
prob_s = ones(1, length(H));
posterior = ones(1, length(H));
for i = 1:1:N
    toss = data(i);
    
    for j = 1:1:length(H)
        posterior(j) = H(j)^toss * (1 - H(j))^(1 - toss) * prior(H(j));
    end
    
    prob_s = prob_s .* posterior;
    prob_s = prob_s / trapz(H, prob_s);
end

% One-step data analysis: easier, but working only for small number of
% tossees (N < 1000)

R = sum(data); % number of heads

% normalization condition ~ int(prob(H)dH) = 1
const = factorial(R) * factorial(N-R) / factorial(N+1);

posterior_1s =@(H) likelihood(R, N, H) .* prior(H) / const;

% calculating probability distrubution

% One-step
prob_1s = ones(1, length(H));
for i = 1:1:length(H)
    prob_1s(i) = posterior_1s(H(i));
end

% check normalization
norm_s = trapz(H, prob_s)
norm_1s = trapz(H, prob_1s)

% Statistics evaluation
[~, index_s] = max(prob_s);
[~, index_1s] = max(prob_1s);

figure()
plot(H, prob_s, 'k')
hold on
plot(H, prob_1s, 'r')
xline(H(index_s), 'k');
xline(H(index_1s), 'r');
xline(w, 'b');

title('Posterior pdf')
xlabel('H')
ylabel('P(H|{data},I)')
xlim([0 1])
legend('pdf_{sequential}', 'pdf_{one-shot}', 'expected value_{seq}', ...
       'expected value_{1s}', 'real bias')


% Prior -------------------------------------------------------------------
function P = prior (H)
    if H >= 0 && H <= 1
        P = 1;
    else
        P = 0;
    end
end

% Likelihood --------------------------------------------------------------
function P = likelihood (R, N, H)
    P = H.^R .* (1 - H).^(N - R);
end