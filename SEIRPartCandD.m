%%% MA 493 - Project 2
% Name: Srinivas Ganesan
% Date: 04/21/2023

%PART 2)c): SEIR Epidemiological Model

% Set the final time 
tFinal = 100;

% Set the step size in the centered difference approximation (5)
h = 0.1/30;

% Create the time vector with interval size dt
t = [0:h:tFinal];

% Set the nominal values of the parameters
k1 = 1e-5;
k2 = 0.03;
k3 = 0.01;
S_init = 1e6; % S_init = S(0)
E_init = 0; % E_init = E(0)
I_init = 0.001*S_init; % I_init = I(0)
R_init = 0; % R_init = R(0)

% 'p' is a list containing the nominal values of the parameters
 p(1) = k1;
 p(2) = k2;
 p(3) = k3;
 p(4) = S_init;
 p(5) = E_init;
 p(6) = I_init;
 p(7) = R_init;

% Set the value of the perturbation for computing scaled sensitivities
alpha = 0.1;

%(II) 
% A) Calculate the scaled sensitivities via the centered difference rule for Susceptibility(S)(6)
% Use the ode solver repeatedly to get (disp1 and disp1f), (disp2 and disp2f), (disp3 and disp3f), and (disp4 and disp4f)
% corresponding to s1,s2,s3,s4
% Solve for the approximated scaled sensitivities using
% Centered finite difference method. Obtain a caliberated 'h' value where
% the sensitivity plot "visually" seems to converge


% (i) Calculate (disp1 and disp1f) for 's1':

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(SEIR(t,y,p+[h,0,0,0,0,0,0]));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);

% Extract the Susceptibility component from 'y'
disp1 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[h,0,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp1f = y(:,1);
s1 = 2*alpha*mNom*(disp1-disp1f)/(2*h);

% (ii) Calculate (disp2 and disp2f) for 's2':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,h,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp2 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,h,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp2f = y(:,1);
s2 = 2*alpha*mNom*(disp2-disp2f)/(2*h);


% (iii) Calculate (disp3 and disp3f) for 's3':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,h,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp3 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,h,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp3f = y(:,1);
s3 = 2*alpha*mNom*(disp3-disp3f)/(2*h);

% (iv) Calculate (disp4 and disp4f) for 's4':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,h,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)+h p(5) p(6) p(7)], options);
disp4 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,h,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)-h p(5) p(6) p(7)], options);
disp4f = y(:,1);
s4 = 2*alpha*mNom*(disp4-disp4f)/(2*h);

% (v) Calculate (disp5 and disp5f) for 's5':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,h,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5)+h p(6) p(7)], options);
disp5 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,h,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5)-h p(6) p(7)], options);
disp5f = y(:,1);
s5 = 2*alpha*mNom*(disp5-disp5f)/(2*h);

% (vi) Calculate (disp6 and disp6f) for 's6':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,0,h,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6)+h p(7)], options);
disp6 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,0,h,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6)-h p(7)], options);
disp6f = y(:,1);
s6 = 2*alpha*mNom*(disp6-disp6f)/(2*h);

% (vii) Calculate (disp7 and disp7f) for 's7':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,0,0,h]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)+h], options);
disp7 = y(:,1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,0,0,h]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)-h], options);
disp7f = y(:,1);
s7 = 2*alpha*mNom*(disp7-disp7f)/(2*h);

% Plot the scaled sensitivities
figure(4)
subplot(2,4,1)
plot(t,s1);
title('dS/dk1')
hold on
subplot(2,4,2);
plot(t,s2);
title('dS/dk2')
hold on
subplot(2,4,3);
plot(t,s3);
title('dS/dk3')
hold on
subplot(2,4,4);
plot(t,s4);
title('dS/dS_init')
hold on
subplot(2,4,5);
plot(t,s5);
title('dS/dE_init')
hold on
subplot(2,4,6);
plot(t,s6);
title('dS/dI_init')
hold on
subplot(2,4,7);
plot(t,s7);
title('dS/dR_init')
hold on

%(III) 
% A) Calculate the scaled sensitivities via the centered difference rule for Exposed count(E)(6)
% Use the ode solver repeatedly to get (disp1 and disp1f), (disp2 and disp2f), (disp3 and disp3f), and (disp4 and disp4f)
% corresponding to s1,s2,s3,s4
% Solve for the approximated scaled sensitivities using
% Centered finite difference method. Obtain a caliberated 'h' value where
% the sensitivity plot "visually" seems to converge


% (i) Calculate (disp1 and disp1f) for 's1':

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(SEIR(t,y,p+[h,0,0,0,0,0,0]));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);

% Extract the 'Exposed' component from 'y'
disp1 = y(:,2);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[h,0,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp1f = y(:,2);
s1 = 2*alpha*mNom*(disp1-disp1f)/(2*h);

% (ii) Calculate (disp2 and disp2f) for 's2':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,h,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp2 = y(:,2);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,h,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp2f = y(:,2);
s2 = 2*alpha*mNom*(disp2-disp2f)/(2*h);


% (iii) Calculate (disp3 and disp3f) for 's3':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,h,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp3 = y(:,2);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,h,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp3f = y(:,2);
s3 = 2*alpha*mNom*(disp3-disp3f)/(2*h);

% (iv) Calculate (disp4 and disp4f) for 's4':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,h,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)+h p(5) p(6) p(7)], options);
disp4 = y(:,2);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,h,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)-h p(5) p(6) p(7)], options);
disp4f = y(:,2);
s4 = 2*alpha*mNom*(disp4-disp4f)/(2*h);

% (v) Calculate (disp5 and disp5f) for 's5':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,h,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5)+h p(6) p(7)], options);
disp5 = y(:,2);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,h,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5)-h p(6) p(7)], options);
disp5f = y(:,2);
s5 = 2*alpha*mNom*(disp5-disp5f)/(2*h);

% (vi) Calculate (disp6 and disp6f) for 's6':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,0,h,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6)+h p(7)], options);
disp6 = y(:,2);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,0,h,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6)-h p(7)], options);
disp6f = y(:,2);
s6 = 2*alpha*mNom*(disp6-disp6f)/(2*h);

% (vii) Calculate (disp7 and disp7f) for 's7':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,0,0,h]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)+h], options);
disp7 = y(:,2);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,0,0,h]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)-h], options);
disp7f = y(:,2);
s7 = 2*alpha*mNom*(disp7-disp7f)/(2*h);

% Plot the scaled sensitivities
figure(5)
subplot(2,4,1)
plot(t,s1);
title('dE/dk1')
hold on
subplot(2,4,2);
plot(t,s2);
title('dE/dk2')
hold on
subplot(2,4,3);
plot(t,s3);
title('dE/dk3')
hold on
subplot(2,4,4);
plot(t,s4);
title('dE/dS_init')
hold on
subplot(2,4,5);
plot(t,s5);
title('dE/dE_init')
hold on
subplot(2,4,6);
plot(t,s6);
title('dE/dI_init')
hold on
subplot(2,4,7);
plot(t,s7);
title('dE/dR_init')
hold on

%(IV) 
% A) Calculate the scaled sensitivities via the centered difference rule for Infected count(I)(6)
% Use the ode solver repeatedly to get (disp1 and disp1f), (disp2 and disp2f), (disp3 and disp3f), and (disp4 and disp4f)
% corresponding to s1,s2,s3,s4
% Solve for the approximated scaled sensitivities using
% Centered finite difference method. Obtain a caliberated 'h' value where
% the sensitivity plot "visually" seems to converge


% (i) Calculate (disp1 and disp1f) for 's1':

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(SEIR(t,y,p+[h,0,0,0,0,0,0]));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);

% Extract the 'Exposed' component from 'y'
disp1 = y(:,3);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[h,0,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp1f = y(:,3);
s1 = 2*alpha*mNom*(disp1-disp1f)/(2*h);

% (ii) Calculate (disp2 and disp2f) for 's2':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,h,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp2 = y(:,3);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,h,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp2f = y(:,3);
s2 = 2*alpha*mNom*(disp2-disp2f)/(2*h);


% (iii) Calculate (disp3 and disp3f) for 's3':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,h,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp3 = y(:,3);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,h,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp3f = y(:,3);
s3 = 2*alpha*mNom*(disp3-disp3f)/(2*h);

% (iv) Calculate (disp4 and disp4f) for 's4':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,h,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)+h p(5) p(6) p(7)], options);
disp4 = y(:,3);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,h,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)-h p(5) p(6) p(7)], options);
disp4f = y(:,3);
s4 = 2*alpha*mNom*(disp4-disp4f)/(2*h);

% (v) Calculate (disp5 and disp5f) for 's5':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,h,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5)+h p(6) p(7)], options);
disp5 = y(:,3);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,h,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5)-h p(6) p(7)], options);
disp5f = y(:,3);
s5 = 2*alpha*mNom*(disp5-disp5f)/(2*h);

% (vi) Calculate (disp6 and disp6f) for 's6':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,0,h,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6)+h p(7)], options);
disp6 = y(:,3);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,0,h,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6)-h p(7)], options);
disp6f = y(:,3);
s6 = 2*alpha*mNom*(disp6-disp6f)/(2*h);

% (vii) Calculate (disp7 and disp7f) for 's7':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,0,0,h]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)+h], options);
disp7 = y(:,3);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,0,0,h]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)-h], options);
disp7f = y(:,3);
s7 = 2*alpha*mNom*(disp7-disp7f)/(2*h);

% Plot the scaled sensitivities
figure(6)
subplot(2,4,1)
plot(t,s1);
title('dI/dk1')
hold on
subplot(2,4,2);
plot(t,s2);
title('dI/dk2')
hold on
subplot(2,4,3);
plot(t,s3);
title('dI/dk3')
hold on
subplot(2,4,4);
plot(t,s4);
title('dI/dS_init')
hold on
subplot(2,4,5);
plot(t,s5);
title('dI/dE_init')
hold on
subplot(2,4,6);
plot(t,s6);
title('dI/dI_init')
hold on
subplot(2,4,7);
plot(t,s7);
title('dI/dR_init')
hold on

%(V) 
% A) Calculate the scaled sensitivities via the centered difference rule for Recovered count(E)(6)
% Use the ode solver repeatedly to get (disp1 and disp1f), (disp2 and disp2f), (disp3 and disp3f), and (disp4 and disp4f)
% corresponding to s1,s2,s3,s4
% Solve for the approximated scaled sensitivities using
% Centered finite difference method. Obtain a caliberated 'h' value where
% the sensitivity plot "visually" seems to converge


% (i) Calculate (disp1 and disp1f) for 's1':

% Set options on for the ode solver
options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 

% Specify the model to be solved - see below for function definition
modelRHS = @(t,y)(SEIR(t,y,p+[h,0,0,0,0,0,0]));

% Solve the ODE system, numerically, using "ode45"
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);

% Extract the 'Exposed' component from 'y'
disp1 = y(:,4);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[h,0,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp1f = y(:,4);
s1 = 2*alpha*mNom*(disp1-disp1f)/(2*h);

% (ii) Calculate (disp2 and disp2f) for 's2':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,h,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp2 = y(:,4);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,h,0,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp2f = y(:,4);
s2 = 2*alpha*mNom*(disp2-disp2f)/(2*h);


% (iii) Calculate (disp3 and disp3f) for 's3':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,h,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp3 = y(:,4);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,h,0,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)], options);
disp3f = y(:,4);
s3 = 2*alpha*mNom*(disp3-disp3f)/(2*h);

% (iv) Calculate (disp4 and disp4f) for 's4':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,h,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)+h p(5) p(6) p(7)], options);
disp4 = y(:,4);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,h,0,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4)-h p(5) p(6) p(7)], options);
disp4f = y(:,4);
s4 = 2*alpha*mNom*(disp4-disp4f)/(2*h);

% (v) Calculate (disp5 and disp5f) for 's5':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,h,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5)+h p(6) p(7)], options);
disp5 = y(:,4);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,h,0,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5)-h p(6) p(7)], options);
disp5f = y(:,4);
s5 = 2*alpha*mNom*(disp5-disp5f)/(2*h);

% (vi) Calculate (disp6 and disp6f) for 's6':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,0,h,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6)+h p(7)], options);
disp6 = y(:,4);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,0,h,0]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6)-h p(7)], options);
disp6f = y(:,4);
s6 = 2*alpha*mNom*(disp6-disp6f)/(2*h);

% (vii) Calculate (disp7 and disp7f) for 's7':

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p+[0,0,0,0,0,0,h]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)+h], options);
disp7 = y(:,4);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-8); 
modelRHS = @(t,y)(SEIR(t,y,p-[0,0,0,0,0,0,h]));
[t,y] = ode45(modelRHS, [0:h:tFinal], [p(4) p(5) p(6) p(7)-h], options);
disp7f = y(:,4);
s7 = 2*alpha*mNom*(disp7-disp7f)/(2*h);

% Plot the scaled sensitivities
figure(7)
subplot(2,4,1)
plot(t,s1);
title('dR/dk1')
hold on
subplot(2,4,2);
plot(t,s2);
title('dR/dk2')
hold on
subplot(2,4,3);
plot(t,s3);
title('dR/dk3')
hold on
subplot(2,4,4);
plot(t,s4);
title('dR/dS_init')
hold on
subplot(2,4,5);
plot(t,s5);
title('dR/dE_init')
hold on
subplot(2,4,6);
plot(t,s6);
title('dR/dI_init')
hold on
subplot(2,4,7);
plot(t,s7);
title('dR/dR_init')
hold on

% Functon "SEIR"
%
% Obtains solution to the SEIR epidemiological model
% Inputs:
%   t - time
%   y - state variables (S(Susceptible),E(Exposed), I(Infected),R(Recovered))
%   p - vector of model parameters [k1, k2, k3, S_init, E_init, I_init, R_init]
% Note: S_init = S(0), E_init = E(0), I_init = I(0), R_init = R(0). They
% are the initial conditions in the SEIR model
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