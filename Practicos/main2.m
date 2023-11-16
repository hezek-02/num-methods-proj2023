##p1  ej 6

h= 0.1.^(1:20);
a=1;
g= @(x) sin(5*x);
gPrima= 5*cos(5*a);
derivadas = zeros(1,20);
for k=1:20
  derivadas(k) = (g(a+h(k))-g(a-h(k)))/(2*h(k));
end
errores = abs(derivadas - gPrima(a)*ones(1,20));
loglog(h,errores,'bo-s','linewidth',1);
xlabel('h');
ylabel('Error absoluto');
title('cocienteIncrCentrado vs derivada, sen(5x)');
grid on;

%%dasd
A=[12,6;3,7];
cond(A);