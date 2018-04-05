% PhD in Economics and Finance
% Bocconi University
% A.y. 2017-2018
%
% Macroeconomics 4
% TA Session 1: The Basic RBC Model
% Andrea Pasqualini

%--------------------------------------------------------------------------
% 1. Defining variables - first endogenous (var) and then exogenous (varexo)
%--------------------------------------------------------------------------

var
    y c n k w r i g a             % endogenous variables
    y_l c_l n_l i_l w_l k_l a_l g_l;  % log variables

varexo
    e_g e_a;                          % exogenous variables


%--------------------------------------------------------------------------
% 2. A. Set the parameters - names
%--------------------------------------------------------------------------

parameters
    beta sigma phi delta alpha
    g_ss rho_g s_g a_ss rho_a s_a;   % stochastic processes


%--------------------------------------------------------------------------
% 2. B. Include the values assigned to the parameters
%--------------------------------------------------------------------------

beta = 0.99;   % discount factor
phi = 1;       % inverse of frisch elasticity
sigma = 1;     % risk aversion (if it is equal to one utility is log)
delta = 0.025; % depreciation rate
alpha = 0.33;  % fraction of output that capital owners receive

rho_a = 0.9;   % persistence of technologic shocks
a_ss = 0.452;  % technology steady state (to have Y=1)
s_a = 0.01;    % standard deviation of technology shocks

rho_g = 0.9;   % persistence of public spending shocks
g_ss = 0.2;    % public spending steady state (to have G/Y=20%)
s_g = 0.01;    % standard deviation of public spending shocks


%--------------------------------------------------------------------------
% 3. Model - Insert the model (if in nonlinear form remember that you will
%            need to include the initial block value (steady state),
%            otherwise you can insert it in linear form,after linearization)
%--------------------------------------------------------------------------

model;

    n^(phi)*c^(sigma) = w;                             % labor supply
    c^(-sigma )= (c(1)^(-sigma))*beta*(1+r(1)-delta);  % euler equation
    y = a*(k(-1)^(alpha))*n^(1-alpha);                 % production function
    r = alpha*y/k(-1);                                 % capital demand
    w = (1-alpha)*y/n;                                 % labor demand
    y = c + i + g;                                     % resource constraint
    k = (1-delta)*k(-1) + i;                           % law of motion of cap
    g_l = (1-rho_g)*log(g_ss)+rho_g*g_l(-1)+e_g;       % gov spending process
    a_l = (1-rho_a)*log(a_ss)+rho_a*a_l(-1)+e_a;       % technology process

    % log variables
    y_l = log(y);    i_l=log(i);      n_l=log(n);
    c_l = log(c);    w_l = log(w);    k_l=log(k);
    a_l=log(a);      g_l =log(g);

end;


%--------------------------------------------------------------------------
% 4.  specify the initial values for the steady state algorithm. You could
% also insert an additional .m file that computes the steady state values
%--------------------------------------------------------------------------

initval;
    y   =  1.00;
    c   =  0.56;
    n   =  1.08;
    k   =  9.43;
    w   =  0.61;
    r   =  0.03;
    g   =  0.20;
    a   =  0.45;
    i   =  0.23;
    y_l =  0.00;
    c_l = -0.56;
    n_l =  0.08;
    i_l = -1.44;
    w_l = -0.48;
    a_l = -0.79;
    g_l = -1.60;
    e_g =  0;
    e_a =  0;
end;


%--------------------------------------------------------------------------
% 5. Compute the Steady State and check if BK conditions are satisfied
%--------------------------------------------------------------------------

steady (solve_algo=0);
check;


%--------------------------------------------------------------------------
% 6.  Now decide which exogenous values you want to shock and by how much
%--------------------------------------------------------------------------

shocks;
    var e_g = s_g^2; % variance of public spending shock
    var e_a = s_a^2; % variance of technologic shock
end;


%--------------------------------------------------------------------------
% 7. Stochastic simulation
%--------------------------------------------------------------------------

stoch_simul (order=1,irf=60)

% variables of interest
y_l
k_l
c_l
n_l
i_l
r
w_l
a_l
g_l
;
