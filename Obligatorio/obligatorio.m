###VARIABLES GLOBALES/CTES###
warning('off', 'all'); %Sacar las warning de los events, ojo saca todas las warnings, pero no teniamos maś
tol = 1e-9
itMax = 250;
c = 299792.458; %km/s




###PARTE 1###
fprintf("\n \n PARTE 1 \n \n");
%Coordenadas de los satelites y sus tiempos de transmisión
          %x        %y         %z           %t
S = [15600,   7540,   20140,   0.07074 ;  %satelite 1
       18760,   2750,   18610,   0.07220  ;  %satelite 2
       17610,   14630, 13480,   0.07690  ;  %satelite 3
       19170,   610,     18390,   0.07242  ];  %satelite 4

F = @(x) sistema_ecs(x,S,0);
JF = @(x) jacobiana(x,S);

x0 = [0, 0, 6370, 0];

[res1,k1] = newton_raphson(x0, JF, F, tol, itMax);
fprintf("\nResultado newton_raphson: (x=%d  y=%d  z=%d  d=%d)\n",res1(1),res1(2),res1(3),res1(4));
fprintf("cant de iteraciones: %d \n \n", k1);
###FIN PARTE1###



###PARTE 2###
fprintf("\n \n PARTE 2 \n \n");

rho = 26570; %km
d = 0.0001; % corrección temporal
fi = [0, pi/6, pi/4, pi/3];
theta = [pi/3, pi/4, 3*pi/4, pi]; 

                        %x                                                 %y                              %z            %t 
S = [ rho*cos(fi(1))*cos(theta(1)),   rho*cos(fi(1))*sin(theta(1)),   rho*sin(fi(1)),    0.0; %satelite 1
         rho*cos(fi(2))*cos(theta(2)),   rho*cos(fi(2))*sin(theta(2)),   rho*sin(fi(2)),   0.0; %satelite 2
         rho*cos(fi(3))*cos(theta(3)),   rho*cos(fi(3))*sin(theta(3)),   rho*sin(fi(3)),   0.0; %satelite 3
         rho*cos(fi(4))*cos(theta(4)),   rho*cos(fi(4))*sin(theta(4)),   rho*sin(fi(4)),   0.0]; %satelite 4

R = zeros(1,4);
t  = zeros (1,4);
for i=1:4
    R(1,i) = sqrt( S(i,1).^2 + S(i,2).^2 + (S(i,3) - 6370).^2 );
    t(1,i) = d + R(1,i) / c;
    S(i,4) = t(1,i);
endfor

F = @(x) sistema_ecs(x,S,d);
JF = @(x) jacobiana(x,S);

x0 = [0, 0, 6370];

[res2,k2] = newton_raphson(x0, JF, F, tol, itMax);
res2(4) = d; % le añado su correccion temporal dado por letra
fprintf("\nResultado newton_raphson: (x=%d  y=%d  z=%d  d=%d) \n",res2(1),res2(2),res2(3),res2(4));
fprintf("cant de iteraciones: %d", k2);
distancia_pto_real = sqrt( (res2(1) - x0(1))^2 + (res2(2) - x0(2))^2 + (res2(3) - x0(3))^2 );
fprintf("\ndistancia_pto_real:%d \n",distancia_pto_real);

delta_ti = perturbaciones(1);
errores_de_salida = zeros(1,16);
factores_de_incremento = zeros(1,16);
for i=1:16
    S(1,4) = t(1,1) + delta_ti(i,1);
    S(2,4) = t(1,2) + delta_ti(i,2);
    S(3,4) = t(1,3) + delta_ti(i,3);
    S(4,4) = t(1,4) + delta_ti(i,4);

    F = @(x) sistema_ecs(x,S,d);
    JF = @(x) jacobiana(x,S);

    x0 = [0, 0, 6370];

    [res3,k3] = newton_raphson(x0, JF, F, tol, itMax);
##    fprintf("\nNewton_raphson perturbado,  con delta_ti=(%d, %d, %d, %d): \n" , delta_ti(i,1), delta_ti(i,2), delta_ti(i,3), delta_ti(i,4));
##
##    fprintf("resultado: (x=%d  y=%d  z=%d)\n", res3(1),res3(2),res3(3));
##    fprintf("cant de iteraciones: %d\n", k3);

    errores_de_salida(1,i) = norm([res2(1) - res3(1),res2(2) - res3(2),res2(3) - res3(3)],inf); %‖(∆x, ∆y, ∆z)‖∞
    factores_de_incremento(1,i) = errores_de_salida(1,i)/(c*10e-8);
endfor

max_err_salida = max(errores_de_salida);
num_cond = max(factores_de_incremento);
fprintf("Dado perturbaciones realizadas a valores de ti en cada satélite hallamos el  mayor error de salida: %d\n ",max_err_salida);
fprintf("y el numero de condicion: %d\n ", num_cond);
###FIN PARTE 2###



###PARTE 3###
fprintf("\n \n PARTE 3 \n \n");

%Coordenadas de los satelites mas juntos
         %x        %y         %z             %t
S = [133,    6326,    25157,   0.0 ;  %satelite 1
       -5269,   7172,   25035,   0.0 ;  %satelite 2
       -5038,   8349,   24716,   0.0 ;  %satelite 3
       -4800,   7698,   24973,   0.0 ];  %satelite 4
       
R_j = zeros(1,4);
t_j  = zeros (1,4);

for i=1:4
    R_j(1,i) = sqrt( S(i,1).^2 + S(i,2).^2 + (S(i,3) - 6370).^2 );
    t_j(1,i) = d + R_j(1,i) / c; 
    S(i,4) = t_j(1,i);
endfor

x0 = [0, 0, 6370];
F = @(x) sistema_ecs(x,S,d);
JF = @(x) jacobiana(x,S);

[res2,k2] = newton_raphson(x0, JF, F, tol, itMax);
res2(4) = d; 

fprintf("\nResultado newton_raphson: (x=%d  y=%d  z=%d  d=%d) \n",res2(1),res2(2),res2(3),res2(4));
fprintf("cant de iteraciones: %d", k2);

distancia_pto_real = sqrt( (res2(1) - x0(1))^2 + (res2(2) - x0(2))^2 + (res2(3) - x0(3))^2 );
fprintf("\ndistancia_pto_real:%d \n",distancia_pto_real);


errores_de_salida = zeros(1,16);
factores_de_incremento = zeros(1,16);
for i=1:16
  S(1,4) = t_j(1,1) + delta_ti(i,1);
  S(2,4) = t_j(1,2) + delta_ti(i,2);
  S(3,4) = t_j(1,3) + delta_ti(i,3);
  S(4,4) = t_j(1,4) + delta_ti(i,4);
  
  F = @(x) sistema_ecs(x,S,d);
  JF = @(x) jacobiana(x,S);

  x0 = [0, 0, 6370];
  
  [res3,k3] = newton_raphson(x0, JF, F, tol, itMax);

##  fprintf("\nNewton_raphson perturbado,  con delta_ti=(%d, %d, %d, %d): \n" , delta_ti(i,1), delta_ti(i,2), delta_ti(i,3), delta_ti(i,4));
##
##  fprintf("resultado: (x=%d  y=%d  z=%d)\n", res3(1),res3(2),res3(3));
##  fprintf("cant de iteraciones: %d\n", k3);

  errores_de_salida(1,i) = norm([res2(1) - res3(1),res2(2) - res3(2),res2(3) - res3(3)],inf); %‖(∆x, ∆y, ∆z)‖∞
  factores_de_incremento(1,i) = errores_de_salida(1,i)/(c*10e-8);
endfor

max_err_salida = max(errores_de_salida);
num_cond = max(factores_de_incremento);

fprintf("\nDado perturbaciones realizadas a valores de ti en cada satélite hallamos el  mayor error de salida: %d\n ",max_err_salida);
fprintf("y el numero de condicion: %d\n ", num_cond);
###FIN PARTE 3 ###



###PARTE 4 ###
fprintf("\n \n PARTE 4 \n \n");

d = 0.0001;

S = [133,           347,        26567;
         24214,     -7393,      8061  ;
         -25087,    -5690,      6651  ;
         -24025,    2219,       11129;
         -10347,    -16476,   18096;
         4948,        19967,    16816;
         -6251,       25465,    4291  ;
         -13410,     -21235,  8673] ;

R = zeros(1,8);
t  = zeros(1,8);
for i =1:8
    R(1,i) = sqrt( S(i,1).^2 + S(i,2).^2 + (S(i,3) - 6370).^2 );
    t(1,i) = d + R(1,i) / c;
endfor

F = @(x) sistema_ecs(x,S,d);
JF = @(x) jacobiana(x,S);
 
%xSol =   [0, 0, 3670, d];

y = zeros(8,1);
y = c * t(1,1:8)';

x0 = [0, 0, 6370, 0];

[res2,k2] = gauss_newton(x0,y, JF, F, tol, itMax);

fprintf("\nResultado gauss_newton: (x=%d  y=%d  z=%d  d=%d) \n",res2(1),res2(2),res2(3),res2(4));
fprintf("cant de iteraciones: %d \n", k2);

distancia_pto_real = sqrt( (res2(1) - x0(1))^2 + (res2(2) - x0(2))^2 + (res2(3) - x0(3))^2 );
fprintf("\ndistancia_pto_real:%d \n",distancia_pto_real);

delta_ti = perturbaciones(2);
errores_de_salida_minCuadrados = zeros(1,256);
factores_de_incremento_minCuadrados = zeros(1,256);
for i=1:256
    y (1,1) = c * (t(1,1) + delta_ti(i,1));
    y (2,1) = c * (t(1,2) + delta_ti(i,2));
    y (3,1) = c * (t(1,3) + delta_ti(i,3));
    y (4,1) = c * (t(1,4) + delta_ti(i,4));
    y (5,1) = c * (t(1,5) + delta_ti(i,5));
    y (6,1) = c * (t(1,6) + delta_ti(i,6));
    y (7,1) = c * (t(1,7) + delta_ti(i,7));
    y (8,1) = c * (t(1,8) + delta_ti(i,8));

    F = @(x) sistema_ecs(x,S,d);
    JF = @(x) jacobiana(x,S);

    [res3,k3] = gauss_newton(x0, y, JF, F, tol, itMax);
    
##    fprintf("\nNewton_raphson perturbado,  con delta_ti=(%d, %d, %d, %d, %d, %d, %d, %d): \n" , delta_ti(i,1), delta_ti(i,2), delta_ti(i,3), delta_ti(i,4),delta_ti(i,5), delta_ti(i,6), delta_ti(i,7), delta_ti(i,8));
##    fprintf("resultado: (x=%f  y=%f  z=%f)\n", res3(1), res3(2), res3(3), res3(4));
##    fprintf("cant de iteraciones: %d\n", k3);

    errores_de_salida_minCuadrados(1,i) = norm([res2(1) - res3(1),res2(2) - res3(2),res2(3) - res3(3), res2(4) - res3(4)],inf); %‖(∆x, ∆y, ∆z)‖∞
    factores_de_incremento_minCuadrados(1,i) = errores_de_salida_minCuadrados(1,i)/(c*10e-8);
endfor
max_err_salida = max(errores_de_salida_minCuadrados);
num_cond = max(factores_de_incremento_minCuadrados);
fprintf("\nDado perturbaciones realizadas a valores de ti en cada satélite hallamos el  mayor error de salida: %d\n ", max_err_salida);
fprintf("y el numero de condicion: %d\n ", num_cond);

###Comparación  de errores perturbados, caso satelites juntos y muchos satélites 8 (min cuad)
% Gráfica de errores perturbados
##figure(1)
##
##plot(1:16, log(errores_de_salida(1:16)+1), '-b.', 1:256, log(errores_de_salida_minCuadrados(1:256)+1), '-r.');
##title('Comparación de Errores Perturbados');
##xlabel('Perturbación');
##ylabel('Log(Error + 1)');
##legend('Errores ', 'Errores Mínimos Cuadrados');
##
##% Gráfica de factores de incremento
##figure(2)
##plot(1:16, log(factores_de_incremento(1:16)+1), '-b.', 1:256, log(factores_de_incremento_minCuadrados(1:256)+1), '-r.');
##title('Comparación de Factores de Increm ento del error');
##xlabel('Perturbación');
##ylabel('Log(Factor de Incremento + 1)');
##legend('Factores ', 'Factores Mínimos Cuadrados');
##
####figure(1)
##plot(1:16, errores_de_salida(1:16), '-b.',1:256, errores_de_salida_minCuadrados(1:256), '-r.');
##figure(2)
##plot( 1:16, factores_de_incremento(1:16), '-b.', 1:256, factores_de_incremento_minCuadrados(1:256),'-r.');
###FIN PARTE 4 ###




## PARTE 5 ##
fprintf("\n \n PARTE 5 \n \n");


d = 0.0001;
c = 299792.458; %km/s
rho = 26570; %km

t_5s = [0.0678669,0.0735332,0.0835988,0.1031279,0.1093036,0.0990338];
t_0s = [0.0678670,0.0735332,0.0835987,0.1031277,0.1093034,0.0990337];
theta = [pi/3,7*pi/11,4*pi/5,9*pi/8,4*pi/3,16*pi/9];

        %x                  %y        
S = [rho*cos(theta(1)), rho*sin(theta(1)); %satelite 1
     rho*cos(theta(2)), rho*sin(theta(2)); %satelite 2
     rho*cos(theta(3)), rho*sin(theta(3)); %satelite 3
     rho*cos(theta(4)), rho*sin(theta(4)); %satelite 4
     rho*cos(theta(5)), rho*sin(theta(5)); %satelite 5
     rho*cos(theta(6)), rho*sin(theta(6))]; %satelite 6

F   = @(x) sistema_ecs(x,S,d);
JF  = @(x) jacobiana(x,S);
y   = zeros(6,1);

x0 = [0,0,0];

y = c * t_5s(1,1:6)';
[res1,k1] = gauss_newton(x0,y, JF, F, tol, itMax);

y = c * t_0s(1,1:6)';
[res2,k2] = gauss_newton(x0,y, JF, F, tol, itMax);

fprintf("\nLa ubicación del globo en -5s es: (x=%.15d  y=%.15d  d=%.15d) \n",res1(1),res1(2),res1(3));
fprintf("cant de iteraciones: %d \n", k1);

fprintf("\nLa ubicación del globo en 0s es: (x=%.15d  y=%.15d  d=%.15d) \n",res2(1),res2(2),res2(3));
fprintf("cant de iteraciones: %d \n", k2);
###FIN PARTE5###



###PARTE 6###
fprintf("\n \n PARTE 6 \n \n");
%valores iniciales
rho_globo_5s = sqrt(res1(1)^2 + res1(2)^2);
rho_globo_0s = sqrt(res2(1)^2 + res2(2)^2);

theta_globo_5s = atan(res1(2)/res1(1));
theta_globo_0s = atan(res2(2)/res2(1));

rho_prima_inicial = (rho_globo_0s - rho_globo_5s)/5;
theta_prima_inicial = (theta_globo_0s - theta_globo_5s)/5;

g = @(rho) 3.986e5./ (rho^2);  %km^3/s^2
k = 1; %cte
M = 1; %1kg;

%y(1)  = rho
%y(2) = theta
%y(3) = rho'
%y(4) = theta'
y0 = [rho_globo_0s, theta_globo_0s, rho_prima_inicial, theta_prima_inicial]';

F_tY_t = @(t,y)  [  
    y(3);  % rho'
    y(4);  % theta'
    -g(y(1)) - (k/M)  *   v(y)  *   y(3) + y(1)  *  y(4).^2      ; %rho''
    -(k/M)    * v(y)    *   y(4)  -   2 * ( y(3) * y(4) )  /  y(1)   ]; %theta''
opciones = odeset('Events', @enSuperficie);
% Resolver la ecuación diferencial con ode23
[t, sol] = ode23(F_tY_t, [0 1000], y0,opciones);

% Graficar resultados
figure(1);
subplot(2, 1, 1);
plot(t, sol(:, 1),  'color', [0, 0.5, 0], 'LineWidth', 2);  % Rho
xlabel('Tiempo en seg');
ylabel('Rho en km');
title('Solución de la Ecuación Diferencial: Rho');

subplot(2, 1, 2);
plot(t, sol(:, 2), 'b', 'LineWidth', 2);  % Theta
xlabel('Tiempo en seg');
ylabel('Theta en radianes');
title('Solución de la Ecuación Diferencial: Theta');

fprintf("\nLa velocidad de impacto del globo con la superficie terrestre es aproximado: %f  km/s   ( %f km/h),\n el tiempo que demoró en colisionar desde t=0 es aproximadamente:  %f  (segundos),  con un ángulo de %f  radianes ( %f  grados)\n",v(sol(end,:)),v(sol(end,:))*3600, t(end),  sol(end,2), sol(end,2)*180/pi);

###FIN PARTE6###



###PARTE 7###
fprintf("\n \n PARTE 7 \n \n");

% se q k se encuentra entre 3.17
kl= 2;
kh=4;
it =0;
velocidad= 0;
tol = 1e-12;
vel_tope = 200;
while (velocidad <vel_tope- tol ||  velocidad >vel_tope- eps)  && it < itMax 
   it = it+1;
   k = (kh +kl)/2; % analogo al metodo de bisección
   
    F_tY_t2 = @(t,y)  [  
        y(3);  % rho'
        y(4);  % theta'
        -g(y(1)) - (k/M)  *   v(y)  *   y(3) + y(1)  *  y(4).^2      ; %rho''
        -(k/M)    * v(y)    *   y(4)  -   2 * ( y(3) * y(4) )  /  y(1)   ]; %theta''
        
    sol =0;
    [t, sol] = ode23(F_tY_t2, [0 800], y0,opciones);
    
    velocidad = v(sol(end,:)) * 3600;
    rho = sol(end,1);
    dif = abs(velocidad - vel_tope-tol);

    if (velocidad> vel_tope && it < 500 )
        kl = k;
    elseif (velocidad< vel_tope- tol && it < 500)
        kh = k;
    endif
endwhile

fprintf("\nLa velocidad de impacto del globo con la superficie terrestre es aproximado: %f  km/s   ( %f km/h),\n el tiempo que demoró en colisionar desde t=0 es aproximadamente:  %f  (segundos),  con un ángulo de %f  radianes ( %f  grados)\n",v(sol(end,:)),v(sol(end,:))*3600, t(end),  sol(end,2), sol(end,2)*180/pi);
fprintf("El valor de k tal que al colisionar con la superficie terrestre el golpe es menor o igual a 200km/h  es: %f \n", k);

###FIN PARTE 7###
