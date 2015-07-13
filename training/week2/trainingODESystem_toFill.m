% Plots the time series of a very simple 6-species network
function trainingODESystem_toFill
% Plots the time series of a very simple 6-species network
% The network includes the following reactions:
%   A + B -> C + D          (with rate constant k1)
%   C -> E + F              (with rate constant k2)
%   D -> B                  (with rate constant k3);
%
% Based on example code by Brian Ingalls:
% http://www.math.uwaterloo.ca/~bingalls/MMSB/Code/matlab/network_example.m

% Rate constant parameters
k1 = NaN; % mM/sec
k2 = NaN; % 1/sec
k3 = NaN; % 1/sec
parameters = [k1 k2 k3];

% Initial species concentrations
S0 = NaN;

% Set ODE simulation parameters
ODE=@trainingODEdt;
options=odeset('Refine', 6);
Tend=10;

% Run simulation
% Use ode45 function and provide it with initial conditions and parameters

% Plot time series
% Fill in

end

function dS = trainingODEdt( ~, s, p )
% Simple ODE system for training purposes

    % Production and degradation terms
    k1 = p(1);
    % Fill in others
    
    % Differential chemical species variables
    a=s(1); 
    % Fill in others

    % Differential System
    a_dt = NaN;                    % fill in;
    b_dt = NaN;                    % fill in;
    c_dt = NaN;                    % fill in;
    d_dt = NaN;                    % fill in;
    e_dt = NaN;                    % fill in;
    f_dt = NaN;                    % fill in;
    
    dS=[a_dt,b_dt,c_dt,d_dt,e_dt,f_dt]';
end