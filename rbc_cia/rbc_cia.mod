% PhD in Economics and Finance
% Bocconi University
% A.y. 2017-2018
%
% Macroeconomics 4
% TA Session 1: The RBC Model with Money (cash in advance constraint)
% Andrea Pasqualini

%--------------------------------------------------------------------------
% 1. Defining variables - first endogenous (var) and then exogenous (varexo)
%--------------------------------------------------------------------------

var
    Y       % output
    A       % labor productivity
    C       % consumption
    N       % labor
    b       % real bond holdings
    m       % real money holdings
    lambda  % budget multiplier
    mu      % cia multiplier
    w       % real wage
    i       % nominal interest rate
    pi      % inflation rate
    r       % real interest rate
    g       % money growth rate
    
    % variables in logs
    Y_l
    A_l
    C_l
    N_l
    b_l
    m_l
    w_l
    i_l
    pi_l
    r_l
    g_l;
    
varexo
    e_a     % productivity shock
    e_g;    % money growth shock

%--------------------------------------------------------------------------
% 2A. Set the parameters - names
%--------------------------------------------------------------------------

parameters
    beta      % discount factor
    gamma     % risk aversion
    vphi      % inverse of Frisch elasticity
    chi       % leisure utility weight
    A_ss      % technology in steady state
    g_ss      % money growth in steady state
    theta     % persistence of technology shocks
    phi       % persistence of money growth shocks
    sigma_a   % standard deviation of productivity shocks
    sigma_g;  % standard deviation of money growth shocks


%--------------------------------------------------------------------------
% 2B. Include the values assigned to the parameters (calibration)
%--------------------------------------------------------------------------

beta = 0.99;   % discount factor
vphi = 1;      % inverse of frisch elasticity
gamma = 1;     % risk aversion (if it is equal to one utility is log)
chi = set_chi; % leisure utility weight

theta = 0.9;     % persistence of technologic shocks
A_ss = 1;        % technology steady state (to have Y=1)
sigma_a = 0.01;  % standard deviation of technology shocks

phi = 0.9;        % persistence of money growth
g_ss = 1.005;     % money growth in steady state
sigma_g = 0.001;  % standard deviation of money growth shocks


%--------------------------------------------------------------------------
% 3. Model - Insert the model (if in nonlinear form remember that you will
%            need to include the initial block value (steady state),
%            otherwise you can insert it in linear form, after linearization)
%--------------------------------------------------------------------------

model;

    w = A;
    Y = A * N;
    C^(-gamma) = lambda + mu;
    chi * N^(vphi) = lambda * w;
    lambda = beta * lambda(+1) * (1 + r(+1));
    lambda = beta * lambda(+1) * (1 / (1 + pi(+1))) + mu;
    C + b + m = w * N + (1 + r) * b(-1) + m(-1) / (1 + pi);
    C = m;
    Y = C;
    1 + r = (1 + i) / (1 + pi);
    m = (g / (1 + pi)) * m(-1);
    log(A) = (1 - theta) * log(A_ss) + theta * log(A(-1)) + e_a;
    log(g) = (1 - phi) * log(g_ss) + phi * log(g(-1)) + e_g;
    
    % variables in logs
    Y_l = log(Y);
    A_l = log(A);
    C_l = log(C);
    N_l = log(N);
    b_l = log(b);
    m_l = log(m);
    w_l = log(w);
    i_l = log(i);
    pi_l = log(pi);
    r_l = log(r);
    g_l = log(g);

end;


%--------------------------------------------------------------------------
% 4. Specify the initial values for the steady state algorithm. You could
% also insert an additional .m file that computes the steady state values
%--------------------------------------------------------------------------

rbc_cia_steady

initval;

    Y      = Y_ss;  % output
    A      = A_ss;  % labor productivity
    C      = C_ss;  % consumption
    N      = N_ss;  % labor
    b      = b_ss;  % real bond holdings
    m      = m_ss;  % real money holdings
    lambda = lambda_ss;  % budget multiplier
    mu     = mu_ss;  % cia multiplier
    w      = w_ss;  % real wage
    i      = i_ss;  % nominal interest rate
    pi     = pi_ss;  % inflation rate
    r      = r_ss;  % real interest rate
    g      = g_ss;  % money growth gross rate
    
    e_g = 0;
    e_a = 0;
    
    % log-variables
    Y_l      = log(Y);
    A_l      = log(A);
    C_l      = log(C);
    N_l      = log(N);
    b_l      = log(b);
    m_l      = log(m);
    w_l      = log(w);
    i_l      = log(1 + i) - 1;
    pi_l     = log(1 + pi) - 1;
    r_l      = log(1 + r) - 1;
    g_l      = log(g);

end;

% resid;

%--------------------------------------------------------------------------
% 5. Compute the Steady State and check if BK conditions are satisfied
%--------------------------------------------------------------------------

steady (solve_algo=0);
check;


%--------------------------------------------------------------------------
% 6. Now decide which exogenous values you want to shock and by how much
%--------------------------------------------------------------------------

shocks;

    var e_a = sigma_a^2;  % variance of technologic shock
    var e_g = sigma_g^2;  % variance of money growth shock

end;


%--------------------------------------------------------------------------
% 7. Stochastic simulation
%--------------------------------------------------------------------------

stoch_simul (order=1, irf=60)
% variables of interest
Y_l
b_l
C_l
N_l
i_l
pi_l
r_l
w_l
m_l
A_l
g_l;

