%funcion a estimar raiz la cual es log3
f =@(x) exp(x) -  3;
f_prima=@(x) exp(x);
%g
g = @(x,d) x - f(x)/d;

x0 = 1;
itMax = 50;
res = ones(3,1);
tol = 1e-12;
x_iter = ones(6,3); %matriz con las iteraciones de cada d 
fprintf('\nDe antemano conociendo el valor de la raíz tendré como criterio de corte |x_k -log(3)| < 1e-9\n y que las iteraciones sean al menos 6, para así poder calcular el cociente del error como se plantea \n. si no, 50 iteraciones máximo\n\n');
for d = 1:3
    x = x0;
    k = 0;
    while k < itMax && (abs(x-log(3)) > tol || k < 6)
      k = k+1;      
      if (k <= 6) 
        x_iter(k,d) = x; %para parte del error
      endif
      x = g(x, d);
    endwhile
    res(d) = x;
    fprintf('resultado de emplear la iteracion con d=%d, nro de iteraciones=%d, convergió a: %f\n\n',d, k, res(d));
endfor

pts = linspace(0,5,150);
figure(1);
 plot(0:5, x_iter(1:6,1), '-ro',0:5, x_iter(1:6,2), '-bo',0:5, x_iter(1:6,3), '-go', pts, log(3), '-k.');
 legend('d=1',"d=2","d=3", "f(log(3)))=0");
 title("primeras 6 iteraciones de g(x,d)=x-f(x)/d");
xlabel('iteraciones(k)');
ylabel('g(x,d) d=1 vs d=2 vs d=3 vs solución');

figure(2);
 plot(0:5, x_iter(1:6,2), '-bo',0:5, x_iter(1:6,3), '-go',pts, log(3),'-k.');
 legend("d=2","d=3", "f(log(3))=0");
 title("primeras 6 iteraciones de g(x,d)=x-f(x)/d");
xlabel('iteraciones(k)');
ylabel('g(x,d) d=2 vs d=3 vs solución');

%ordenes estimados
p = [1, 2];

disp("ordenes estimados;  para d=2 es 1 y orden para d=3 es 2, (las que convergen) ");

error_k = zeros (2,5);
for d =2:3
    fprintf("\nCaso d = %d \n", d);
    x = x0;
    for k=1:5 
        ek = abs(x_iter(k+1,d) -  log(3));
        ekp = abs(x_iter(k,d) -  log(3))^p(d-1);
        error_k(d-1,k) = ek/ekp;
        fprintf('error en cociente con k = %d, %f  \n', k-1, error_k(d-1,k));
    endfor    
endfor

g_prima = @(x,d) 1 - f_prima(x)/d;
g_prima_prima = @(x,d) - f_prima(x)/d;


err_d2= g_prima(res(2),1)/1;% p=1
fprintf('\nvelocidad de convergencia d=2 y p=1:  %f \n', err_d2);

##err_d2 = g_prima_prima(res(2),2)/2; % p=2
##fprintf('velocidad de convergencia d=2 y p=2 %f \n', err_d2);

err_d3= g_prima(res(3),3)/1; % p=1
fprintf('velocidad de convergencia d=3 y p=1:  %f \n', err_d3);

err_d3 = g_prima_prima(res(3),3)/2; %p=2
fprintf('velocidad de convergencia d=3 y p=2:  %f \n', err_d3);


figure(3);
 plot(0:4, log(error_k(1, 1:5)+1), '-bo', 0:4,log(error_k(2, 1:5)+1) , '-go');
 legend('Err coc.  d=2','Err. coc. d=3');
 title('Error cociente  d=2 vs d=3 (en escala logarítmica)');
xlabel('k');
ylabel('Error');

## Guardar la primera figura
##figure(1);
##print('figure1.png', '-dpng', '-r300'); 
##
#### Guardar la segunda figura
##figure(2);
##print('figure2.png', '-dpng', '-r300');
##
#### Guardar la tercera figura
##figure(3);
##print('figure3.png', '-dpng', '-r300');

## figure;
##plot(0:4, error_k(1, 1:5), '-bo');
##title('Error para d=2');
##xlabel('k');
##ylabel('Error');
##
##figure;
##plot(0:4, error_k(2, 1:5), '-go');
##title('Error para d=3');
##xlabel('k');
##ylabel('Error');


