%% MA 493 - Mathematical Foundations of Data Science
%% Response function for Topic IV.2 - equation (4) in the notes

function f = fOD(t,theta)

m = theta(1);
k = theta(2);
c = theta(3);
A = theta(4);

D1 = sqrt(c^2 - 4*k*m);
f = A/(2*D1)*((c + D1)*exp(-(c - D1)/(2*m)*t) - (c - D1)*exp(-(c + D1)/(2*m)*t));