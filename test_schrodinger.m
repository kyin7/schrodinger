function [psi, opt] = test_schrodinger
%% domain parameters
opt.k = 0.0001; % time step size
opt.T = 1; % maximum time
opt.t = 0:opt.k:opt.T;
opt.x1 = -8:1/32:8; % one space dimension
opt.x2 = -8:1/32:8; % the other space dimension
[opt.X2, opt.X1] = meshgrid(opt.x1, opt.x2);
opt.n1 = length(opt.x1);
opt.n2 = length(opt.x2);

%% parameter setting: O(1)-interactions, positive initial energy
% opt.epsilon = 1.0;
% opt.kappa = -2.0;
% opt.gamma1 = 1.0;
% opt.gamma2 = 1.0;
% opt.V = reshape((opt.gamma1^2*opt.X1.^2 + opt.gamma2^2*opt.X2.^2)/2, [], 1); % potential
% opt.psi0 = 1/sqrt(pi*opt.epsilon) * exp(-opt.V/opt.epsilon); % initial state
% opt.mu = 1e2; % L1 penalty parameter

%% Strong interactions, negative initial energy
% opt.epsilon = 0.3;
% opt.kappa = -1.9718;
% opt.gamma1 = 1.0;
% opt.gamma2 = 1.0;
% V = (opt.gamma1^2*opt.X1.^2 + opt.gamma2^2*opt.X2.^2)/2; % potential
% opt.V = reshape(V, [], 1);
% psi0 = 1/sqrt(pi*opt.epsilon) * exp(-opt.V/opt.epsilon); % initial state
% opt.psi0 = reshape(psi0, [], 1);
% opt.mu = 1e2; % L1 penalty parameter

%% Very weak interactions, anisotropic condensate, nonzero initial phase
% opt.epsilon = 1.0;
% opt.kappa = 0.1;
% opt.gamma1 = 1.0;
% opt.gamma2 = 2.0;
% V = (opt.gamma1^2*opt.X1.^2 + opt.gamma2^2*opt.X2.^2)/2; % potential
% opt.V = reshape(V, [], 1);
% psi0 = opt.gamma2^0.25/sqrt(pi*opt.epsilon) * exp(-V/opt.epsilon) ...
%     .* exp(1j*cosh(sqrt(opt.X1.^2+2*opt.X2.^2))/opt.epsilon); % initial state
% opt.psi0 = reshape(psi0, [], 1);
% opt.mu = 1e2; % L1 penalty parameter

%% Strong interactions, nonzero initial phase data
opt.epsilon = 0.1;
opt.kappa = 1.2649;
opt.gamma1 = 1.0;
opt.gamma2 = 1.0;
opt.mus = sqrt(opt.kappa*opt.gamma2/pi);
V = (opt.gamma1^2*opt.X1.^2 + opt.gamma2^2*opt.X2.^2)/2; % potential
opt.V = reshape(V, [], 1);
psi0 = (opt.X1.^2+opt.X2.^2<2*opt.mus)...
    .*sqrt((opt.mus-(opt.X1.^2+opt.X2.^2)/2)/opt.kappa)...
    .*exp(1j/opt.epsilon*cosh(sqrt(opt.X1.^2+2*opt.X2.^2))); % initial state
opt.psi0 = reshape(psi0, [], 1);
opt.mu = 1; % L1 penalty parameter

%% run the algorithm
[psi] = solve_schrodinger_L1(opt);
