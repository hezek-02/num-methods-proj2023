
function [detener,isterm,dir] = enSuperficie(t,y)
    detener = y(1) - 6370;  % Detiene la integración cuando rho =6370
    isterm = 1; % que termine si se cumple la condicion
    dir = -1; % que la condicion se alcance desde abajo