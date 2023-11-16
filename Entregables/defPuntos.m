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
## @deftypefn {} {@var{retval} =} defPuntos (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ezequiel Capó <hezek@heze-Makina>
## Created: 2023-10-01

function [x,y] = defPuntos(grado)
  %modo chebyschev o normal
  f = @(x) 1./(1 + 25*x.^2);
  n = grado+1; %qty puntos interp 
  x=zeros(1,n);
  %evaluo abscisa equidistante (-1,1) y sus imagenes
  for j=1:n
    x(j)= -1 + 2*(j-1)./grado; %equivalente a linespace
  endfor  
  %evalua valores f(x)
  for i=1:n
    y(i)= f(x(i));
  end
endfunction
