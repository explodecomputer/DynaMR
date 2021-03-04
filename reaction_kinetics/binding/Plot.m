clear all
clc

global a b Kt

a = 1;
b = 0.5;

% Time interval (min)
ts = 0.1;
t = 0:ts:30;

% Options;
options = optimset('TolFun',1e-12);

% Initial values
Xt = 1;
Kt = 0.3;

Y0_full = [Xt Kt 0];

% Solving the Full system
[T,Full] = ode23s(@Sys_Full,t,Y0_full,options);

% Solving the MM system
% [T,MM] = ode23s(@Sys_MM,t,Y0_MM,options);
figure(1)
plot(T,Full(:,1),T,Full(:,2),T,Full(:,3),'linewidth',1)
xlabel('Time','fontsize',18)
ylabel('Pop','fontsize',18)
legend('X','K','XK','location','northeast');
saveas(gcf,'Comparison','png');

