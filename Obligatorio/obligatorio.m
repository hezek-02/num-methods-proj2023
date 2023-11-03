
###PARTE 1###
tic
%Coordenadas de los satelites y sus tiempos de transmisión
          %x        %y         %z           %t
S = [15600,   7540,   20140,   0.07074  ];  %satelite 1
       [18760,   2750,   18610,   0.07220  ];  %satelite 2
       [17610,   14630, 13480,   0.07690  ];  %satelite 3
       [19170,   610,     18390,   0.07242  ];  %satelite 4

%cte vel de la luz 
c = 299792.458; %km/s

%Sistema de ecs. para hallar la posición del receptor y su correción temporal ( tiempo que demora en recibir la señal)
F = @(x, y, z, d) [
    sqrt((x - S(1,1))).^2  +  (y - S(1,2)).^2    + (z - S(1,3)).^2)   -   c * (S(1,4) - d) ;
    sqrt((x - S(2,1))).^2  +  (y - S(2,2))).^2   + (z - S(2,3)).^2)   -   c * (S(2,4) - d) ;
    sqrt((x - S(3,1))).^2  +  (y - S(3,2))).^2   + (z - S(3,3)).^2)   -   c * (S(3,4) - d) ;
    sqrt((x - S(4,1))).^2  +  (y - S(4,2))).^2   + (z - S(4,3)).^2)   -   c * (S(4,4) - d) ;
];

%jacobina de F, como ejemplo de momento, cambiar por func jacobiana(x,y,z,d) que retorne el resultado, (la calcule a mano)
Jf = @(x, y, z, d) [ 
  (x-S(1,1))/sqrt((x - S(1,1)).^2 + (y - S(1,2)).^2 + (z - S(1,3)).^2),      (y-S(1,2))/sqrt((x - S(1,1)).^2 + (y - S(1,2)).^2 + (z - S(1,3)).^2),        (z-S(1,3))/sqrt((x - S(1,1)).^2 + (y - S(1,2)).^2 + (z - S(1,3)).^2),       c ;                        
  (x-S(2,1))/sqrt((x - S(2,1)).^2 + (y - S(2,2)).^2 + (z - S(2,3)).^2),      (y-S(2,2))/sqrt((x - S(2,1)).^2 + (y - S(2,2)).^2 + (z - S(2,3)).^2),        (z-S(2,3))/sqrt((x - S(2,1)).^2 + (y - S(2,2)).^2 + (z - S(2,3)).^2),       c ;
  (x-S(3,1))/sqrt((x - S(3,1)).^2 + (y - S(3,2)).^2 + (z - S(3,3)).^2),      (y-S(3,2))/sqrt((x - S(3,1)).^2 + (y - S(3,2)).^2 + (z - S(3,3)).^2),        (z-S(3,3))/sqrt((x - S(3,1)).^2 + (y - S(3,2)).^2 + (z - S(3,3)).^2),       c ;
  (x-S(4,1))/sqrt((x - S(4,1)).^2 + (y - S(4,2)).^2 + (z - S(4,3)).^2),      (y-S(4,2))/sqrt((x - S(4,1)).^2 + (y - S(4,2)).^2 + (z - S(4,3)).^2),        (z-S(4,3))/sqrt((x - S(4,1)).^2 + (y - S(4,2)).^2 + (z - S(4,3)).^2),       c ;
];
%punto inicial para iteración NR
x0 = [0, 0, 6370, 0];

%obtenemos la matriz Jacobiana
%res = jacobian([A(1);A(2);A(3);A(4)]);

%empleamos NR
[res,k] = newton_raphson(x0, Jf, F, 1000);
fprintf("\n resultado newton_raphson: ");
disp(res);
fprintf("cant de iteraciones: %d", k);
fprintf("\n \n");


###PARTE 2###

%Ai = ρ cos(fi ) cos(theta)
%Bi = ρ cos(fi ) sen(theta)
%Ci = ρ sen(fi )
%0 ≤ fi ≤ π/2 y 0 ≤ theta ≤ 2π
%ϕ1 = 0, ϕ2 = π/6, ϕ3 = π/4, ϕ4 = π/3
rho = 26570 %km, magnitud

fi = [0, pi/6, pi/4, pi/3];
theta = [0, pi/4, 3*pi/4, pi]; %valores distintos dos a dos

%Coordenadas de los satelites y sus tiempos de transmisión
                        %x                                                 %y                              %z                  %t 
S = [  rho*cos(fi(1))*cos(theta(1)),   rho*cos(fi(1))*sen(theta(1)),   rho*sen(fi(1)),   0.091239]; %satelite 1
       [  rho*cos(fi(2))*cos(theta(2)),   rho*cos(fi(2))*sen(theta(2)),   rho*sen(fi(2)),   0.091239]; %satelite 2
       [  rho*cos(fi(3))*cos(theta(3)),   rho*cos(fi(3))*sen(theta(3)),   rho*sen(fi(3)),   0.091239]; %satelite 3
       [  rho*cos(fi(4))*cos(theta(4)),   rho*cos(fi(4))*sen(theta(4)),   rho*sen(fi(4)),   0.091239]; %satelite 4

R = zeros(4,1);
for i=1:4
    R(i) = sqrt( S(i,1)^2 + S(i,2)^2 + (S(i,3) - 6370)^2)

endfor

%jacobiana de F
Jf = @(x, y, z) [
  (x-S1(1))/sqrt((x - S1(1)).^2 + (y - S1(2)).^2 + (z - S1(3)).^2),  (y-S1(2))/sqrt((x - S1(1)).^2 + (y - S1(2)).^2 + (z - S1(3)).^2), (z-S1(3))/sqrt((x - S1(1)).^2 + (y - S1(2)).^2 + (z - S1(3)).^2) ;
  (x-S2(1))/sqrt((x - S2(1)).^2 + (y - S2(2)).^2 + (z - S2(3)).^2),  (y-S2(2))/sqrt((x - S2(1)).^2 + (y - S2(2)).^2 + (z - S2(3)).^2), (z-S2(3))/sqrt((x - S2(1)).^2 + (y - S2(2)).^2 + (z - S2(3)).^2) ;
  (x-S3(1))/sqrt((x - S3(1)).^2 + (y - S3(2)).^2 + (z - S3(3)).^2),  (y-S3(2))/sqrt((x - S3(1)).^2 + (y - S3(2)).^2 + (z - S3(3)).^2), (z-S3(3))/sqrt((x - S3(1)).^2 + (y - S3(2)).^2 + (z - S3(3)).^2) ;
  (x-S4(1))/sqrt((x - S4(1)).^2 + (y - S4(2)).^2 + (z - S4(3)).^2),  (y-S4(2))/sqrt((x - S4(1)).^2 + (y - S4(2)).^2 + (z - S4(3)).^2), (z-S4(3))/sqrt((x - S4(1)).^2 + (y - S4(2)).^2 + (z - S4(3)).^2) ;
];


%punto inicial para iteración NR
x0 = [0, 0, 6370, 0];


%empleamos NR
[res,k] = newton_raphson(x0, Jf, F, 1000);
fprintf("\n resultado newton_raphson: ");
disp(res);
fprintf("cant de iteraciones: %d", k);
fprintf("\n \n");

toc