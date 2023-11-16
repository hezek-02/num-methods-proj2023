%% Ejemplo de uso del solver ode45

clear 
clc
close all

a = 1; % parametro para f
f = @(t,y) a*y*(1-y); % y'(t)= f(t,y(t))
intervalo = [0, 8]; % en que intervalo queremos resolver el PVI
y0 = 0.1; % condicion inicial

% ode45(f, intervalo, y0)

[t,y] = ode45(f, intervalo, y0);
% plot(t,y, '-o')

%% Solvers distintos generan soluciones distintas

close all

ode23s(f, intervalo, y0);
hold on
ode23(f, intervalo, y0);