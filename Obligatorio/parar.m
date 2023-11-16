function [val,isterm,dir] = parar(t,y)
val = y(1)-100; % condicion que define al evento
isterm = 1; % que termine si se cumple la condicion
dir = 1; % que la condicion se alcance desde abajo