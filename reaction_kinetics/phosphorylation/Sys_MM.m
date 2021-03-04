function dydt = Sys_MM(t,y)

global a b c e f g Kt Pt

dydt = zeros(2,1);

X = y(1);
Y = y(2);

XK = Kt*X/((b+c)/a+X);
YP = Pt*Y/((f+g)/e+Y);
K = Kt - XK;
P = Pt - YP;

dydt(1) = -a*X*K + b*XK + g*YP;
dydt(2) = c*XK - e*Y*P + f*YP;

end
