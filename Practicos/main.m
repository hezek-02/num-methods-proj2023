## Copyright (C) 2023 Ezequiel Capó
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} main (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ezequiel Capó <hezek@heze-Makina>
## Created: 2023-08-01

tic
#ej1a
vector = [1 2 3 4 5 6];
normaVec = norma(vector);
matriz = [1 2 3 4 5 6; 3 2 1 3 2 1];
disp("norma euclidea personal: ");
disp(normaVec);
disp("norma euclidea octave: ");
disp(norm(vector));
#raizProx = aproximarRaiz(6,25);
#disp(raizProx);
##ej1b
disp(sustitucion(matriz,vector,3));
disp("  ")
##ej1c
M = [1 2 3 4 5 6; 3 2 1 3 2 1];
M = intercambio(M,1,2);
disp(M);

##ej2
disp(signo(3));
disp(signo(0));
disp(signo(-5));
disp(esPar(3));
disp(esPar(2));
disp(max(2,10));
disp(max(10,2));
disp(" ")
##ej 3
disp(maxModif(M));

##ej 4
e = 1e-6;
disp(serie2(@Sgeometrica,e));
##ej5
#a=fibonacci(25);
#disp(a);
b=fiboAureo(255);
disp(" ")
disp(b);
#ej5
#semillaMASIVa
maxSem=1;
largoMaxActual = length(collatz(maxSem,200));
for i=1:999
  if ( largoMaxActual< length(collatz(i,200)))
    maxSem = i;  
    largoMaxActual = length(collatz(i,200));
  endif
endfor
  #maxSem = 871
  disp(collatz(maxSem,200));
toc




