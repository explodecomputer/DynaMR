clear all
clc

global a b c e f g Kt Pt

a = 7;
b = 8;
c = 1;
e = 4;
f = 6;
g = 1;


% Time interval (min)
ts = 0.1;
t = 0:ts:30;

% Options;
options = optimset('TolFun',1e-12);

% Initial values
Xt = 1;
Kt = 0.3;
Pt = 0.5;

Y0_full = [Xt 0 Kt 0 Pt 0];
Y0_MM = [Xt 0];

% Solving the Full system
[T,Full] = ode23s(@Sys_Full,t,Y0_full,options);

% Solving the MM system
% [T,MM] = ode23s(@Sys_MM,t,Y0_MM,options);
figure(1)
plot(T,Full(:,1),T,Full(:,2),T,Full(:,3),T,Full(:,4),T,Full(:,5),T,Full(:,6),'linewidth',1)
xlabel('Time','fontsize',18)
ylabel('Pop','fontsize',18)
legend('X','Y','K','XK','P','YP','location','northeast');
saveas(gcf,'Comparison','png');

