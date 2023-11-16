##p1  ej 6

%vector con valores de h
h= 0.1.^(1:20);
a=1; %pto a evaluar
g= @(x) sin(5*x); %funcion a evaluar
gPrima= 5*cos(5*a);%derivada de funcion a evaluar
derivadas = zeros(1,20); 
for k=1:20
  derivadas(k) = (g(a+h(k))-g(a-h(k)))/(2*h(k));%calculo cociente incremental centrado
end
errores = abs(derivadas - gPrima(a)*ones(1,20)); %calculo de errores
loglog(h,errores,'b-*','linewidth',1);
xlabel('h');
ylabel('Error absoluto');
title('cociente Incremental Centrado vs derivada, sen(5x)');

% Guardar los datos
%datos = [h', errores']; 
%save('datos_resultados.txt', 'datos', '-ascii');
