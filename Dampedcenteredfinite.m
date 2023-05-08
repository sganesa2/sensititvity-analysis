%%% MA 493 - Project 2
% Name: Srinivas Ganesan
% Date: 04/21/2023

%PART 1: Damped Oscillations
% Example 2)(b)

% Set the final time 
tFinal = 8;

% Set the step size in the centered difference approximation (5)
h = 0.1/5;

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



%(II) 
% Calculate the scaled sensitivities via the centered difference rule (5)
% Use the ode solver repeatedly to get (disp1 and disp1f), (disp2 and disp2f), (disp3 and disp3f), and (disp4 and disp4f)
% corresponding to s1,s2,s3,s4
% Solve for the approximated scaled sensitivities using
% Centered finite difference method. Obtain a caliberated 'h' value where
% the sensitivity plot "visually" seems to converge


% (i) Calculate (disp1 and disp1f) for 's1':

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(Damped(t,y,p+[h,0,0,0]));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) 0], options);

% Extract the displacement component from 'y'
disp1 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(Damped(t,y,p-[h,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) 0], options);
disp1f = y(:,1);
s1 = 2*alpha*mNom*(disp1-disp1f)/(2*h);

% (ii) Calculate ('disp2' and 'disp2f') for 's2':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(Damped(t,y,p+[0,h,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) 0], options);
disp2 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(Damped(t,y,p-[0,h,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) 0], options);
disp2f = y(:,1);
s2 = 2*alpha*mNom*(disp2-disp2f)/(2*h);

% (iii) Calculate ('disp3' and 'disp3f') for 's3':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(Damped(t,y,p+[0,0,h,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) 0], options);
disp3 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(Damped(t,y,p-[0,0,h,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) 0], options);
disp3f = y(:,1);
s3 = 2*alpha*mNom*(disp3-disp3f)/(2*h);

% (iv) Calculate ('disp4' and 'disp4f') for 's4':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(Damped(t,y,p+[0,0,0,h]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)+h 0], options);
disp4 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(Damped(t,y,p-[0,0,0,h]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)-h 0], options);
disp4f = y(:,1);
s4 = 2*alpha*mNom*(disp4-disp4f)/(2*h);


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

% Compute the senstivity matrix based on centered differences

% Calculate (dispData1 and dispData1f), (dispData2 and dispData2f), (dispData3 and dispData3f), and (dispData4 and dispData4f which are
% displacement of the mass at different time given by tData.
% These displacement values correspond to (S:,1), (S:,2), (S:,3), (S:,4)

% (i) Calculate ('dispData1' and 'dispData1f') for '(S,:1)':
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(tData,y)(Damped(tData,y,p+[h,0,0,0]));
[tData,y] = ode45(modelRHS, [2:2:8], [p(4) 0], options);
dispData1 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(tData,y)(Damped(tData,y,p-[h,0,0,0]));
[tData,y] = ode45(modelRHS, [2:2:8], [p(4) 0], options);
dispData1f = y(:,1);

S(:,1) = 2*alpha*mNom*(dispData1-dispData1f)/(2*h);


% (ii) Calculate ('dispData2' and 'dispData2f') for 'S(:,2)':
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(tData,y)(Damped(tData,y,p+[0,h,0,0]));
[tData,y] = ode45(modelRHS, [2:2:8], [p(4) 0], options);
dispData2 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(tData,y)(Damped(tData,y,p-[0,h,0,0]));
[tData,y] = ode45(modelRHS, [2:2:8], [p(4) 0], options);
dispData2f = y(:,1);

S(:,2) = 2*alpha*mNom*(dispData2-dispData2f)/(2*h);

% (iii) Calculate ('dispData3' and 'dispData3f') for 'S(:,3)':
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(tData,y)(Damped(tData,y,p+[0,0,h,0]));
[tData,y] = ode45(modelRHS, [2:2:8], [p(4) 0], options);
dispData3 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(tData,y)(Damped(tData,y,p-[0,0,h,0]));
[tData,y] = ode45(modelRHS, [2:2:8], [p(4) 0], options);
dispData3f = y(:,1);

S(:,3) = 2*alpha*mNom*(dispData3-dispData3f)/(2*h);


% (iv) Calculate ('dispData4' and 'dispData4f') for 'S(:,4)':
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(tData,y)(Damped(tData,y,p+[0,0,0,h]));
[tData,y] = ode45(modelRHS, [2:2:8], [p(4)+h 0], options);
dispData4 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(tData,y)(Damped(tData,y,p-[0,0,0,h]));
[tData,y] = ode45(modelRHS, [2:2:8], [p(4)-h 0], options);
dispData4f = y(:,1);

S(:,4) = 2*alpha*mNom*(dispData4-dispData4f)/(2*h);

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