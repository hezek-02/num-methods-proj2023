%% Modelo de Lotka-Volterra para dinamica predador-presa

clear
clc
close all

alpha = 0.02;
y0 = [100; 10];
F = @(t,y) [2*y(1) - alpha*y(1)*y(2); -y(2)+alpha*y(1)*y(2)];
[t, y] = ode45(F,[0 20],y0);

subplot(2,1,1)
plot(t, y,'-o', 'LineWidth',2)
xlabel('t', 'FontSize',14)
ylabel('poblaciones', 'FontSize',14)
legend('conejos','zorros')

subplot(2,1,2)
plot(y(:,1),y(:,2))
grid on
xlabel('conejos', 'FontSize',14)
ylabel('zorros', 'FontSize',14)

%% Ejemplo de evento: detectar el periodo de la solucion

clear
clc
close all

alpha = 0.02;
y0 = [100 10];
F = @(t,y) [2*y(1) - alpha*y(1)*y(2); -y(2)+alpha*y(1)*y(2)];

opts = odeset('events', @parar);
[t,y,tfinal] = ode45(F,[0 14],y0,opts);

plot(t, y,'-o', 'LineWidth',2)

t(end)
y(end,:)

plot(y(:,1),y(:,2),'-o')
