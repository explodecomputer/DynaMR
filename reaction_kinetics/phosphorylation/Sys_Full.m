function dydt = Sys_Full(t,y)

global a b c e f g

dydt = zeros(6,1);

X = y(1);
Y = y(2);
K = y(3);
XK = y(4);
P = y(5);
YP = y(6);

dydt(1) = -a*X*K + b*XK + g*YP;
dydt(2) = c*XK - e*Y*P + f*YP;
dydt(3) = -a*X*K + (b+c)*XK;
dydt(4) = a*X*K - (b+c)*XK;
dydt(5) = -e*Y*P + (f+g)*YP;
dydt(6) = e*Y*P - (f+g)*YP;

end
