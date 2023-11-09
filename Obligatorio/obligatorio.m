###PARTE 1###
fprintf("\n \n PARTE 1 \n \n");

%Coordenadas de los satelites y sus tiempos de transmisión
          %x        %y         %z           %t
S = [15600,   7540,   20140,   0.07074 ;  %satelite 1
       18760,   2750,   18610,   0.07220  ;  %satelite 2
       17610,   14630, 13480,   0.07690  ;  %satelite 3
       19170,   610,     18390,   0.07242  ];  %satelite 4

%cte vel de la luz 
c = 299792.458; %km/s

%Sistema de ecs. para hallar la posición del receptor y su correción temporal ( tiempo que demora en recibir la señal),x(1)=x,x(2)=y,x(3)=z, x(4) =d
F = @(x) [
    sqrt((x(1) - S(1,1)).^2  +  (x(2) - S(1,2)).^2   + (x(3) - S(1,3)).^2)   -   c * (S(1,4) - x(4)) ;
    sqrt((x(1) - S(2,1)).^2  +  (x(2) - S(2,2)).^2   + (x(3) - S(2,3)).^2)   -   c * (S(2,4) - x(4)) ;
    sqrt((x(1) - S(3,1)).^2  +  (x(2) - S(3,2)).^2   + (x(3) - S(3,3)).^2)   -   c * (S(3,4) - x(4)) ;
    sqrt((x(1) - S(4,1)).^2  +  (x(2) - S(4,2)).^2   + (x(3) - S(4,3)).^2)   -   c * (S(4,4) - x(4)) ;
];

%jacobina de F  4x4
Jf = @(x) [ 
  (x(1)-S(1,1))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),      (x(2)-S(1,2))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),        (x(3)-S(1,3))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),   c ;                        
  (x(1)-S(2,1))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),      (x(2)-S(2,2))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),        (x(3)-S(2,3))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),   c ;
  (x(1)-S(3,1))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),      (x(2)-S(3,2))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),        (x(3)-S(3,3))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),   c ;
  (x(1)-S(4,1))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),      (x(2)-S(4,2))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),        (x(3)-S(4,3))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),   c ;
];

%vector inicial para iteración NR
x0 = [0, 0, 6370, 0];

%empleamos NR
[res1,k1] = newton_raphson(x0, Jf, F, 1000, 1);
fprintf("\nResultado newton_raphson: x=%d  y=%d  z=%d  t=%d \n",res1(1),res1(2),res1(3),res1(4));
%disp(res1);
fprintf("cant de iteraciones: %d \n \n", k1);


###PARTE 2###
fprintf("\n \n PARTE 2 \n \n");

rho = 26570; %km, magnitud

d = 0.0001; %cte de letra: corrección temporal

fi = [0, pi/6, pi/4, pi/3];

theta = [0, pi/4, 3*pi/4, pi]; %valores distintos dos a dos

%Coordenadas de los satelites y sus tiempos de transmisión
                        %x                                                 %y                              %z            %t 
S = [ rho*cos(fi(1))*cos(theta(1)),   rho*cos(fi(1))*sin(theta(1)),   rho*sin(fi(1)),    0.0; %satelite 1
         rho*cos(fi(2))*cos(theta(2)),   rho*cos(fi(2))*sin(theta(2)),   rho*sin(fi(2)),   0.0; %satelite 2
         rho*cos(fi(3))*cos(theta(3)),   rho*cos(fi(3))*sin(theta(3)),   rho*sin(fi(3)),   0.0; %satelite 3
         rho*cos(fi(4))*cos(theta(4)),   rho*cos(fi(4))*sin(theta(4)),   rho*sin(fi(4)),   0.0]; %satelite 4

R = zeros(1,4);
t  = zeros (1,4);
for i=1:4
    R(1,i) = sqrt( S(i,1).^2 + S(i,2).^2 + (S(i,3) - 6370).^2 );
    t(1,i) = d + R(1,i) / c; %tiempo de transmisión de cada satelite
    S(i,4) = t(1,i);
endfor

%Sistema de ecs. para hallar la posición del receptor
F2 = @(x) [
    sqrt((x(1) - S(1,1)).^2  +  (x(2) - S(1,2)).^2   + (x(3) - S(1,3)).^2)   -   c * (S(1,4) - d) ;
    sqrt((x(1) - S(2,1)).^2  +  (x(2) - S(2,2)).^2   + (x(3) - S(2,3)).^2)   -   c * (S(2,4) - d) ;
    sqrt((x(1) - S(3,1)).^2  +  (x(2) - S(3,2)).^2   + (x(3) - S(3,3)).^2)   -   c * (S(3,4) - d) ;
    sqrt((x(1) - S(4,1)).^2  +  (x(2) - S(4,2)).^2   + (x(3) - S(4,3)).^2)   -   c * (S(4,4) - d) ;
];

%jacobina de F2, 4x3
Jf2 = @(x) [ 
  (x(1)-S(1,1))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),      (x(2)-S(1,2))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),        (x(3)-S(1,3))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2);                    
  (x(1)-S(2,1))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),      (x(2)-S(2,2))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),        (x(3)-S(2,3))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2);
  (x(1)-S(3,1))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),      (x(2)-S(3,2))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),        (x(3)-S(3,3))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2);
  (x(1)-S(4,1))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),      (x(2)-S(4,2))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),        (x(3)-S(4,3))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2);
];

%punto inicial para iteración NR
x0 = [0, 0, 6370];

%empleamos NR
[res2,k2] = newton_raphson(x0, Jf2, F2, 1000,0);
res2(4) = d; % le añado su correccion temporal dado por letra
fprintf("\nResultado newton_raphson: x=%d  y=%d  z=%d  t=%d \n",res2(1),res2(2),res2(3),res2(4));
%disp(res2);
fprintf("cant de iteraciones: %d", k2);
fprintf("\n \n");

%posibles permutaciones
delta_ti = [ 10e-8,    10e-8,    10e-8,     10e-8;

                 -10e-8,   10e-8,    10e-8,      10e-8;
                  10e-8,    -10e-8,   10e-8,      10e-8;
                  10e-8,    10e-8,    -10e-8,     10e-8;
                  10e-8,    10e-8,     10e-8,    -10e-8;

                  -10e-8,   -10e-8,   10e-8,      10e-8;
                  -10e-8,   10e-8,    -10e-8,     10e-8;
                  -10e-8,   10e-8,    10e-8,     -10e-8;
                  10e-8,    -10e-8,   -10e-8,     10e-8;
                  10e-8,    -10e-8,    10e-8,    -10e-8;
                  10e-8,    10e-8,    -10e-8,    -10e-8;

                  -10e-8,   -10e-8,   -10e-8,    10e-8;
                  -10e-8,   10e-8,    -10e-8,    -10e-8;
                  -10e-8,   -10e-8,   10e-8,     -10e-8;
                  10e-8,    -10e-8,   -10e-8,    -10e-8;
                  
                  -10e-8,   -10e-8,   -10e-8,    -10e-8];

%t(1,i); sin alterar
error_salida = zeros(16);
factor_de_incremento = zeros(16);
for i=1:16
  S(1,4) = t(1,1) + delta_ti(i,1);
  S(2,4) = t(1,2) + delta_ti(i,2);
  S(3,4) = t(1,3) + delta_ti(i,3);;
  S(4,4) = t(1,4) + delta_ti(i,4);

  F2 = @(x) [
      sqrt((x(1) - S(1,1)).^2  +  (x(2) - S(1,2)).^2   + (x(3) - S(1,3)).^2)   -   c * (S(1,4) - d) ;
      sqrt((x(1) - S(2,1)).^2  +  (x(2) - S(2,2)).^2   + (x(3) - S(2,3)).^2)   -   c * (S(2,4) - d) ;
      sqrt((x(1) - S(3,1)).^2  +  (x(2) - S(3,2)).^2   + (x(3) - S(3,3)).^2)   -   c * (S(3,4) - d) ;
      sqrt((x(1) - S(4,1)).^2  +  (x(2) - S(4,2)).^2   + (x(3) - S(4,3)).^2)   -   c * (S(4,4) - d) ;
  ];

  Jf2 = @(x) [
    (x(1)-S(1,1))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),      (x(2)-S(1,2))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),        (x(3)-S(1,3))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2);
    (x(1)-S(2,1))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),      (x(2)-S(2,2))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),        (x(3)-S(2,3))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2);
    (x(1)-S(3,1))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),      (x(2)-S(3,2))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),        (x(3)-S(3,3))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2);
    (x(1)-S(4,1))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),      (x(2)-S(4,2))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),        (x(3)-S(4,3))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2);
  ];

  [res3,k3] = newton_raphson(x0, Jf2, F2, 1000,1);
fprintf("\nNewton_raphson perturdado,  con delta_ti= %d, %d, %d, %d): \n" , delta_ti(i,1), delta_ti(i,2), delta_ti(i,3), delta_ti(i,4));

fprintf("resultado, x=%d  y=%d  z=%d:\n", res3(1),res3(2),res3(3));
##  disp(res3);
##  fprintf("cant de iteraciones: %d", k3);
##  fprintf("\n \n");

  errores_de_salida(i) = norm([res2(1)-res3(1),res2(2)-res3(2),res2(3)-res3(3)],inf); %‖(∆x, ∆y, ∆z)‖∞
  factores_de_incremento(i) = errores_de_salida(i)/(c*10e-8);

  %luego de calcular cada valor por separado


endfor
max_err_salida = max(errores_de_salida);
num_cond = max(factores_de_incremento);
fprintf("\n mayor error de salida: %d\n ", max_err_salida);
fprintf("numero de condicion: %d\n ", num_cond);

##
#####PARTE 3###
##%Coordenadas de los satelites mas juntos y sus tiempos de transmisión
##          %x        %y         %z           %t
##S_juntos = [133,   6326,   25157,   0.0 ;  %satelite 1
##       -5269,   7172,   25035,   0.0 ;  %satelite 2
##       -5038,   8349, 24716,   0.0 ;  %satelite 3
##       -4800,   7698,     24973,   0.0 ];  %satelite 4
##       
##R_j = zeros(1,4);
##t_j  = zeros (1,4);
##for i=1:4
##    R_j(1,i) = sqrt( S_juntos(i,1).^2 + S_juntos(i,2).^2 + (S_juntos(i,3) - 6370).^2 );
##    t_j(1,i) = d + R_j(1,i) / c; %tiempo de transmisión de cada satelite
##    S_juntos(i,4) = t(1,i);
##endfor
##
##%Sistema de ecs. para hallar la posición del receptor
##%Sistema de ecs. para hallar la posición del receptor y su correción temporal ( tiempo que demora en recibir la señal),x(1)=x,x(2)=y,x(3)=z, x(4) =d
##F2_juntos = @(x) [
##    sqrt((x(1) - S_juntos(1,1)).^2  +  (x(2) - S_juntos(1,2)).^2   + (x(3) - S_juntos(1,3)).^2)   -   c * (S_juntos(1,4) - d) ;
##    sqrt((x(1) - S_juntos(2,1)).^2  +  (x(2) - S_juntos(2,2)).^2   + (x(3) - S_juntos(2,3)).^2)   -   c * (S_juntos(2,4) - d) ;
##    sqrt((x(1) - S_juntos(3,1)).^2  +  (x(2) - S_juntos(3,2)).^2   + (x(3) - S_juntos(3,3)).^2)   -   c * (S_juntos(3,4) - d) ;
##    sqrt((x(1) - S_juntos(4,1)).^2  +  (x(2) - S_juntos(4,2)).^2   + (x(3) - S_juntos(4,3)).^2)   -   c * (S_juntos(4,4) - d) ;
##];
##
##%jacobina de F2, 4x3
##Jf2_juntos = @(x) [ 
##  (x(1)-S_juntos(1,1))/sqrt((x(1) - S_juntos(1,1)).^2 + (x(2) - S_juntos(1,2)).^2 + (x(3) - S_juntos(1,3)).^2),      (x(2)-S_juntos(1,2))/sqrt((x(1) - S_juntos(1,1)).^2 + (x(2) - S_juntos(1,2)).^2 + (x(3) - S_juntos(1,3)).^2),        (x(3)-S_juntos(1,3))/sqrt((x(1) - S_juntos(1,1)).^2 + (x(2) - S_juntos(1,2)).^2 + (x(3) - S_juntos(1,3)).^2);                    
##  (x(1)-S_juntos(2,1))/sqrt((x(1) - S_juntos(2,1)).^2 + (x(2) - S_juntos(2,2)).^2 + (x(3) - S_juntos(2,3)).^2),      (x(2)-S_juntos(2,2))/sqrt((x(1) - S_juntos(2,1)).^2 + (x(2) - S_juntos(2,2)).^2 + (x(3) - S_juntos(2,3)).^2),        (x(3)-S_juntos(2,3))/sqrt((x(1) - S_juntos(2,1)).^2 + (x(2) - S_juntos(2,2)).^2 + (x(3) - S_juntos(2,3)).^2);
##  (x(1)-S_juntos(3,1))/sqrt((x(1) - S_juntos(3,1)).^2 + (x(2) - S_juntos(3,2)).^2 + (x(3) - S_juntos(3,3)).^2),      (x(2)-S_juntos(3,2))/sqrt((x(1) - S_juntos(3,1)).^2 + (x(2) - S_juntos(3,2)).^2 + (x(3) - S_juntos(3,3)).^2),        (x(3)-S_juntos(3,3))/sqrt((x(1) - S_juntos(3,1)).^2 + (x(2) - S_juntos(3,2)).^2 + (x(3) - S_juntos(3,3)).^2);
##  (x(1)-S_juntos(4,1))/sqrt((x(1) - S_juntos(4,1)).^2 + (x(2) - S_juntos(4,2)).^2 + (x(3) - S_juntos(4,3)).^2),      (x(2)-S_juntos(4,2))/sqrt((x(1) - S_juntos(4,1)).^2 + (x(2) - S_juntos(4,2)).^2 + (x(3) - S_juntos(4,3)).^2),        (x(3)-S_juntos(4,3))/sqrt((x(1) - S_juntos(4,1)).^2 + (x(2) - S_juntos(4,2)).^2 + (x(3) - S_juntos(4,3)).^2);
##];
##
##%punto inicial para iteración NR
##x0 = [0, 0, 6370];
##
##
##%empleamos NR
##[res2_juntos,k2_juntos] = newton_raphson(x0, Jf2_juntos, F2_juntos, 1000,0);
##res2_juntos(4) = d;
##fprintf("\n resultado newton_raphson: ");
##disp(res2_juntos);
##fprintf("cant de iteraciones: %d", k2_juntos);
##fprintf("\n \n");
##
##%posibles permutaciones
##S_juntos(1,4) = S_juntos(1,4)(+/-)10e-8;
##S_juntos(2,4) = S_juntos(2,4)(+/-)10e-8;
##S_juntos(3,4) = S_juntos(3,4)(+/-)10e-8;
##S_juntos(4,4) = S_juntos(4,4)(+/-)10e-8;
##
##F2_juntos = @(x) [
##    sqrt((x(1) - S_juntos(1,1)).^2  +  (x(2) - S_juntos(1,2)).^2   + (x(3) - S_juntos(1,3)).^2)   -   c * (S_juntos(1,4) - d) ;
##    sqrt((x(1) - S_juntos(2,1)).^2  +  (x(2) - S_juntos(2,2)).^2   + (x(3) - S_juntos(2,3)).^2)   -   c * (S_juntos(2,4) - d) ;
##    sqrt((x(1) - S_juntos(3,1)).^2  +  (x(2) - S_juntos(3,2)).^2   + (x(3) - S_juntos(3,3)).^2)   -   c * (S_juntos(3,4) - d) ;
##    sqrt((x(1) - S_juntos(4,1)).^2  +  (x(2) - S_juntos(4,2)).^2   + (x(3) - S_juntos(4,3)).^2)   -   c * (S_juntos(4,4) - d) ;
##];
##
##Jf2_juntos = @(x) [
##  (x(1)-S_juntos(1,1))/sqrt((x(1) - S_juntos(1,1)).^2 + (x(2) - S_juntos(1,2)).^2 + (x(3) - S_juntos(1,3)).^2),      (x(2)-S_juntos(1,2))/sqrt((x(1) - S_juntos(1,1)).^2 + (x(2) - S_juntos(1,2)).^2 + (x(3) - S_juntos(1,3)).^2),        (x(3)-S_juntos(1,3))/sqrt((x(1) - S_juntos(1,1)).^2 + (x(2) - S_juntos(1,2)).^2 + (x(3) - S_juntos(1,3)).^2);
##  (x(1)-S_juntos(2,1))/sqrt((x(1) - S_juntos(2,1)).^2 + (x(2) - S_juntos(2,2)).^2 + (x(3) - S_juntos(2,3)).^2),      (x(2)-S_juntos(2,2))/sqrt((x(1) - S_juntos(2,1)).^2 + (x(2) - S_juntos(2,2)).^2 + (x(3) - S_juntos(2,3)).^2),        (x(3)-S_juntos(2,3))/sqrt((x(1) - S_juntos(2,1)).^2 + (x(2) - S_juntos(2,2)).^2 + (x(3) - S_juntos(2,3)).^2);
##  (x(1)-S_juntos(3,1))/sqrt((x(1) - S_juntos(3,1)).^2 + (x(2) - S_juntos(3,2)).^2 + (x(3) - S_juntos(3,3)).^2),      (x(2)-S_juntos(3,2))/sqrt((x(1) - S_juntos(3,1)).^2 + (x(2) - S_juntos(3,2)).^2 + (x(3) - S_juntos(3,3)).^2),        (x(3)-S_juntos(3,3))/sqrt((x(1) - S_juntos(3,1)).^2 + (x(2) - S_juntos(3,2)).^2 + (x(3) - S_juntos(3,3)).^2);
##  (x(1)-S_juntos(4,1))/sqrt((x(1) - S_juntos(4,1)).^2 + (x(2) - S_juntos(4,2)).^2 + (x(3) - S_juntos(4,3)).^2),      (x(2)-S_juntos(4,2))/sqrt((x(1) - S_juntos(4,1)).^2 + (x(2) - S_juntos(4,2)).^2 + (x(3) - S_juntos(4,3)).^2),        (x(3)-S_juntos(4,3))/sqrt((x(1) - S_juntos(4,1)).^2 + (x(2) - S_juntos(4,2)).^2 + (x(3) - S_juntos(4,3)).^2);
##];
##
##[res3_juntos,k3_juntos] = newton_raphson(x0, Jf2_juntos, F2_juntos, 1000,0);
##fprintf("\n resultado newton_raphson perturbado: ");
##disp(res3_juntos);
##fprintf("cant de iteraciones: %d", k3_juntos);
##fprintf("\n \n");
##
##error_salida = norm([res2_juntos(1)-res3_juntos(1),res2_juntos(2)-res3_juntos(2),res2_juntos(3)-res3_juntos(3)],inf)
##factor_de_incremento = error_salida/(c*10e-8)
##
##%luego de calcular cada valor por separado
##errores_de_salida = [0.028021, 0.056069, 0.9209, 0.9244, 0.9175, 0.8962, 0.1628, 0.2225];
##max_err_salida = max(errores_de_salida)
##factores_de_incremento = [0.9347, 1.8703, 30.718, 31.428, 30.603, 29.893, 5.4301, 7.4204];
##num_cond = max(factores_de_incremento)
##       
##       
