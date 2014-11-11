function [psi] = solve_schrodinger_L1(opt)
% solve the schrodinger equation with subgradient of L1 term
% psi_t(x,t) = i*epsilon/2*Laplace psi(x,t) - i/(eps*mu)*p(psi(x,t))
% -i/epsilon*(V(x) + kappa1*|psi(x,t)|^2)*psi(x,t)
% where p(u)=u/|u| if u nonzero; p(u)=0 if u=0.
% assume rectangular domain with uniform grid
% Delta x = h, Delta t = k
k = opt.k; % time step size
N = length(opt.t); % total number of time steps
x1 = opt.x1; % one space dimension
x2 = opt.x2; % the other space dimension
n1 = length(x1);
n2 = length(x2);
epsilon = opt.epsilon;
kappa = opt.kappa;
V = opt.V; % potential
mu = opt.mu; % L1 penalty parameter
psi = opt.psi0;
for n=2:N
    % Strang splitting
    % step 1: solve nonlinear ODE
    % psi_t(x,t) = -i/epsilon*(V(x) + kappa1*|psi(x,t)|^2)*psi(x,t)
    % for time step k/2
    psi_temp1 = exp(-1j/epsilon*(V + kappa*abs(psi).^2)*k/2).*psi;
    % step 2: solve "diffusion" equation
    % psi_t(x,t) = i*epsilon/2*Laplace psi(x,t) - i/(eps*mu)*p(psi(x,t))
    % for time step k, intial data the result in step 1
    % use backward in time implicit scheme
    % (I - i*epsilon*k/2*Laplace)psi(x, t_{n+1}) = psi(x, t_n) -
    % i*k/(epsilon*mu)*p(psi(x, t_n))
    rhs = psi_temp1 - 1j*k/(epsilon*mu)*(psi_temp1./(abs(psi_temp1)+eps));
    psi_temp2 = solve_Laplace(rhs, 1, -1j*epsilon*k/2, n1, n2);
    % step 3: repeat step 1 with intial data the result in step 2
    psi = exp(-1j/epsilon*(V + kappa*abs(psi_temp2).^2)*k/2).*psi_temp2;
end
    
