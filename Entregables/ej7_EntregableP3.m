% Función de Runge
f = @(x) 1./(1 + 25*x.^2);

% Definición de los puntos de interpolación

qtyPts=1000;
pts= linspace(-1, 1, qtyPts);%intervalo de pts a graficar
%Funcion interpolante

%grado 4
%sin chebyschev
[x,y] = defPuntos(4);
interp =  interpolacion_Lagrange(x,y,pts);
interpEnPts =  interpolacion_Lagrange(x,y,x);
%con chebyschev
[x_cheb,y_cheb] = defPuntosChebyschev(4);
interp_cheb =  interpolacion_Lagrange(x_cheb,y_cheb,pts);
interpEnPts_cheb =  interpolacion_Lagrange(x_cheb,y_cheb,x_cheb);
figure(1)
plot(pts, f(pts), 'k-',pts, interp, 'b-',x, interpEnPts, 'b*',pts, interp_cheb, 'r-',x_cheb, interpEnPts_cheb, 'r*');
legend('Función de Runge(f)', 'Interpolación de Runge (interp)', 'interp(x)', 'Interpolación de Runge c/Chebyschev(interpCheb)', 'interpCheb(x)', 'Location', 'South');
title('interpolacion Runge grado 4');

%grado 8
%sin chebyschev
[x,y] = defPuntos(8);
interp =  interpolacion_Lagrange(x,y,pts);
interpEnPts =  interpolacion_Lagrange(x,y,x);
%con chebyschev
[x_cheb,y_cheb] = defPuntosChebyschev(8);
interp_cheb =  interpolacion_Lagrange(x_cheb,y_cheb,pts);
interpEnPts_cheb =  interpolacion_Lagrange(x_cheb,y_cheb,x_cheb);
figure(2)
plot(pts, f(pts), 'k-',pts, interp, 'b-',x, interpEnPts, 'b*',pts, interp_cheb, 'r-',x_cheb, interpEnPts_cheb, 'r*');
legend('Función de Runge(f)', 'Interpolación de Runge (interp)', 'interp(x)', 'Interpolación de Runge c/Chebyschev(interpCheb)', 'interpCheb(x)', 'Location', 'South');
title('interpolacion Runge grado 8');

%grado 12
%sin chebyschev
[x,y] = defPuntos(12);
interp =  interpolacion_Lagrange(x,y,pts);
interpEnPts =  interpolacion_Lagrange(x,y,x);
%con chebyschev
[x_cheb,y_cheb] = defPuntosChebyschev(12);
interp_cheb =  interpolacion_Lagrange(x_cheb,y_cheb,pts);
interpEnPts_cheb =  interpolacion_Lagrange(x_cheb,y_cheb,x_cheb);
figure(3)
plot(pts, f(pts), 'k-',pts, interp, 'b-',x, interpEnPts, 'b*',pts, interp_cheb, 'r-',x_cheb, interpEnPts_cheb, 'r*');
legend('Función de Runge(f)', 'Interpolación de Runge (interp)', 'interp(x)', 'Interpolación de Runge c/Chebyschev(interpCheb)', 'interpCheb(x)', 'Location', 'South');
title('interpolacion Runge grado 12');

%grado 16
%sin chebyschev
[x,y] = defPuntos(16);
interp =  interpolacion_Lagrange(x,y,pts);
interpEnPts =  interpolacion_Lagrange(x,y,x);
%con chebyschev
[x_cheb,y_cheb] = defPuntosChebyschev(16);
interp_cheb =  interpolacion_Lagrange(x_cheb,y_cheb,pts);
interpEnPts_cheb =  interpolacion_Lagrange(x_cheb,y_cheb,x_cheb);
figure(4)
plot(pts, f(pts), 'k-',pts, interp, 'b-',x, interpEnPts, 'b*',pts, interp_cheb, 'r-',x_cheb, interpEnPts_cheb, 'r*');
legend('Función de Runge(f)', 'Interpolación de Runge (interp)', 'interp(x)', 'Interpolación de Runge c/Chebyschev(interpCheb)', 'interpCheb(x)', 'Location', 'South');
title('interpolacion Runge grado 16');

%c
[x,y] = defPuntos(6);
interp = interpolacion_lineal_trozos(x, y, pts)
interpEnPts =  interpolacion_lineal_trozos(x,y,x);
figure(5)
plot(pts, f(pts), 'k-',pts, interp, 'b-',x, interpEnPts, 'b*',pts, interp_cheb, 'r-',x_cheb, interpEnPts_cheb, 'r*');
legend('Función de Runge(f)', 'Interpolación de Runge a trozos (interp)', 'interp(x)', 'Interpolación de Runge c/Chebyschev(interpCheb)', 'interpCheb(x)', 'Location', 'South');
title('interpolación lineal a trozos vs interpolación con nodos Chebyschev');




