%%% MA 493 - Project 2
% Name: Srinivas Ganesan
% Date: 04/21/2023

%PART 2)b): SEIR Epidemiological Model

% Main routine to run the SEIR Epidemiological Model
% 
% Set the nominal parameter values
  k1 = 1e-5;
  k2 = 0.03;
  k3 = 0.01;
  S_init = 1e6; % S_init = S(0)
  E_init = 0; % E_init = E(0)
  I_init = 0.001*S_init; % I_init = I(0)
  R_init = 0; % R_init = R(0)

% Specify the final time and the time interval and spacing
tFinal = 100;
dt = 1;
t = [0:dt:tFinal];

% Use ode45 to solve our system of ODEs where:
% t is the vector of time values 
% y is the vector of solutions at the time points in t

% Specify the paramater values (passed to the ODE solver)
  p(1) = k1;
  p(2) = k2;
  p(3) = k3;
  p(4) = S_init;
  p(5) = E_init;
  p(6) = I_init;
  p(7) = R_init;

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(SEIR(t,y,p));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:dt:tFinal], [p(4) p(5) p(6) p(7)], options);

% Extract the relevant solutions
S = y(:,1);
E = y(:,2);
I = y(:,3);
R = y(:,4);

% Plot the solutions
figure(1)
subplot(2,2,1)
plot(t,S,'LineWidth', 3)
title('Susceptible')
hold on

subplot(2,2,2)
plot(t,E,'LineWidth', 3)
title('Exposed')
hold on

subplot(2,2,3)
plot(t,I,'LineWidth', 3)
title('Infected')
hold on

subplot(2,2,4)
plot(t,R,'LineWidth', 3)
title('Recovered')
hold on


% Functon "SEIR"
%
% Obtains solution to the SEIR epidemiological model
% Inputs:
%   t - time
%   y - state variables (S(Susceptible),E(Exposed), I(Infected),R(Recovered))
%   p - vector of model parameters [k1, k2, k3, S_init, E_init, I_init, R_init]
% Output:
%   dy/dt 

function dy = SEIR(t, y, p)
  dy = [0; 0; 0; 0];
  k1 = p(1);
  k2 = p(2);
  k3 = p(3);
  S_init = p(4);
  E_init = p(5);
  I_init = p(6);
  R_init = p(7);

  dy(1) = -k1*y(1)*y(3);
  dy(2) = (k1*y(1)*y(3))-(k2*y(2)) ;
  dy(3) = (k2*y(2))-(k3*y(3));
  dy(4) = k3*y(3);

end