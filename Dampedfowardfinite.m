%%% MA 493 - Project 2
% Name: Srinivas Ganesan
% Date: 04/21/2023

%PART 1: Damped Oscillations
% Example 2)(a)

% Set the final time 
tFinal = 8;

% Set the step size in the forward difference approximation (5)
h = 0.1/20;

% Create the time vector with interval size h
t = [0:h:tFinal];

% Set the nominal values of the parameters
mNom = 1;
kNom = 1;
cNom = 2.3;
ANom = 0.1;

% 'p' is a list containing the nominal values of the parameters
p(1) = mNom;
p(2) = cNom;
p(3) = kNom;
p(4) = ANom;

% Set the value of the perturbation for computing scaled sensitivities
alpha = 0.1;

%(I) Solve the system of 1st order ODE's to get the values for our
% Response function 'y' at different Time's given by 't'

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(Damped(t,y,p));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) 0], options);

% Extract the displacement component from 'y'
dispMain = y(:,1);


%(II) Using the displacement given by dispMain, Solve for the approximated scaled sensitivities using
% Forward finite difference method. Obtain a caliberated 'h' value where
% the sensitivity plot "visually" seems to converge

% Calculate the scaled sensitivities via the forward difference rule (5)
% Use the ode solver repeatedly to get disp1, disp2, disp3, and disp4
% corresponding to s1,s2,s3,s4

% (i) Calculate 'disp1' for 's1':

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(Damped(t,y,p+[h,0,0,0]));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) 0], options);

% Extract the displacement component from 'y'
disp1 = y(:,1);
s1 = 2*alpha*mNom*(disp1-dispMain)/h;

% (ii) Calculate 'disp2' for 's2':

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(Damped(t,y,p+[0,h,0,0]));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) 0], options);

% Extract the displacement component from 'y'
disp2 = y(:,1);
s2 = 2*alpha*mNom*(disp2-dispMain)/h;

% (iii) Calculate 'disp3' for 's3':

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(Damped(t,y,p+[0,0,h,0]));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) 0], options);

% Extract the displacement component from 'y'
disp3 = y(:,1);
s3 = 2*alpha*mNom*(disp3-dispMain)/h;

% (iv) Calculate 'disp4' for 's4':

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(Damped(t,y,p+[0,0,0,h]));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)+h 0], options);

% Extract the displacement component from 'y'
disp4 = y(:,1);
s4 = 2*alpha*mNom*(disp4-dispMain)/h;


% Plot the scaled sensitivities
figure(2)
subplot(2,2,1)
plot(t,s1);
title('dy/dm')
hold on
subplot(2,2,2);
plot(t,s2);
title('dy/dk')
hold on
subplot(2,2,3);
plot(t,s3);
title('dy/dc')
hold on
subplot(2,2,4);
plot(t,s4);
title('dy/dA')
hold on

% (III) Lastly, compute the scaled sensitivity matrix at the 4 indicated 
% data points

% Here, we have 4 data points (rows) and 4 parameters (columns)
S = zeros(4,4);

% Store the time values where data is collected
tData = [2 4 6 8];

% Compute the senstivity matrix based on forward differences

% Calculate dispData, dispData1, dispData2, dispData3, dispData4 which are
% displacement of the mass at different time given by tData.
% These displacement values correspond to yData, (S:,1), (S:,2), (S:,3), (S:,4)

% Calculate dispData
% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(tData,y)(Damped(tData,y,p));

% Solve the ODE system, numerically, using "ode45"
[tData,y] = ode45(modelRHS, [2:2:8], [p(4) 0], options);

% Extract the displacement component from 'y'
dispData = y(:,1);

% a) Calculate dispData1
% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(tData,y)(Damped(tData,y,p+[h,0,0,0]));

% Solve the ODE system, numerically, using "ode45"
[tData,y] = ode45(modelRHS, [2:2:8], [p(4) 0], options);

% Extract the displacement component from 'y'
dispData1 = y(:,1);
S(:,1) = 2*alpha*mNom*(dispData1-dispData)/h;

% b) Calculate dispData2
% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(tData,y)(Damped(tData,y,p+[0,h,0,0]));

% Solve the ODE system, numerically, using "ode45"
[tData,y] = ode45(modelRHS, [2:2:8], [p(4) 0], options);

% Extract the displacement component from 'y'
dispData2 = y(:,1);
S(:,2) = 2*alpha*mNom*(dispData2-dispData)/h;

% c) Calculate dispData3
% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(tData,y)(Damped(tData,y,p+[0,0,h,0]));

% Solve the ODE system, numerically, using "ode45"
[tData,y] = ode45(modelRHS, [2:2:8], [p(4) 0], options);

% Extract the displacement component from 'y'
dispData3 = y(:,1);
S(:,3) = 2*alpha*mNom*(dispData3-dispData)/h;

% d) Calculate dispData4
% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(tData,y)(Damped(tData,y,p+[0,0,0,h]));

% Solve the ODE system, numerically, using "ode45"
[tData,y] = ode45(modelRHS, [2:2:8], [p(4)+h 0], options);

% Extract the displacement component from 'y'
dispData4 = y(:,1);
S(:,4) = 2*alpha*mNom*(dispData4-dispData)/h;

S

% Functon "Damped"
%
% Implements ODE model for damped oscillatory motion of a
% mass-spring-dashpot system

% Inputs:
%   t - time
%   y - state variables (displacement and velocity)
%   p - vector of model parameters [m c k]
% Output:
%   dy/dt 

function dy = Damped(t, y, p)
  dy = [0; 0];
  m = p(1);
  c = p(2);
  k = p(3);

  dy(1) = y(2);
  dy(2) = -k/m * y(1) - c/m * y(2);

end
