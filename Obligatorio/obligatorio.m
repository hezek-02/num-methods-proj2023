
###PARTE 1###
tic
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
%punto inicial para iteración NR
x0 = [0, 0, 6370, 0];

%obtenemos la matriz Jacobiana
%res = jacobian([A(1);A(2);A(3);A(4)]);

%empleamos NR
[res1,k1] = newton_raphson(x0, Jf, F, 1000, 0);
fprintf("\n resultado newton_raphson: ");
disp(res1);
fprintf("cant de iteraciones: %d", k1);
fprintf("\n \n");


###PARTE 2###
fprintf("\n \n PARTE 2 \n \n");

%Ai = ρ cos(fi ) cos(theta)
%Bi = ρ cos(fi ) sen(theta)
%Ci = ρ sen(fi )
%0 ≤ fi ≤ π/2 y 0 ≤ theta ≤ 2π
%ϕ1 = 0, ϕ2 = π/6, ϕ3 = π/4, ϕ4 = π/3

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
%Sistema de ecs. para hallar la posición del receptor y su correción temporal ( tiempo que demora en recibir la señal),x(1)=x,x(2)=y,x(3)=z, x(4) =d
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
res2(4) = d;
fprintf("\n resultado newton_raphson: ");
disp(res2);
fprintf("cant de iteraciones: %d", k2);
fprintf("\n \n");

toc