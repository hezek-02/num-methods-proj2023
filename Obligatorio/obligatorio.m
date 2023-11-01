
###PARTE 1###

%Coordenadas de los satelites y sus tiempos de transmisión
S1 = [15600, 7540, 20140, 0.07074];
S2 = [18760, 2750, 18610, 0.07220];
S3 = [17610, 14630, 13480, 0.07690];
S4 = [19170, 610, 18390, 0.07242];

%cte vel de la luz en 
c = 299792.458; %km/s

%Sistema para hallar la posición del receptor y su correción temporal ( tiempo que demora en recibir la señal)
A = @(x, y, z, d) [
    sqrt((x - S1(1)).^2 + (y - S1(2)).^2 + (z - S1(3)).^2) - c * (S1(4) - d);
    sqrt((x - S2(1)).^2 + (y - S2(2)).^2 + (z - S2(3)).^2) - c * (S2(4) - d);
    sqrt((x - S3(1)).^2 + (y - S3(2)).^2 + (z - S3(3)).^2) - c * (S3(4) - d);
    sqrt((x - S4(1)).^2 + (y - S4(2)).^2 + (z - S4(3)).^2) - c * (S4(4) - d);
];

%punto inicial para iteración NR
x0 = [0, 0, 6370, 0];

%obtenemos la matriz Jacobiana
res = jacobian([A(1);A(2);A(3);A(4)]);

%empleamos NR
res = newton_raphson(x0, )

