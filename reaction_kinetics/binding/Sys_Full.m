function dydt = Sys_Full(t,y)

global a b

dydt = zeros(3,1);

X = y(1);
K = y(2);
XK = y(3);

dydt(1) = -a*X*K + b*XK;
dydt(2) = -a*X*K + b*XK;
dydt(3) = a*X*K - b*XK;

end
