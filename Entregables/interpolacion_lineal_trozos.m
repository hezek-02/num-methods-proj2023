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
## @deftypefn {} {@var{retval} =} interpolacion_linealTrozos (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ezequiel Capó <hezek@heze-Makina>
## Created: 2023-10-01

function w = interpolacion_lineal_trozos(x, y, v)
    % Encontrar el índice i tal que x(i) <= v < x(i+1)
    n=length(x);%cant pts evaluar
    largoV=length(v);
    w=zeros(size(v));
    for i=1:largoV
      for j = 1:n-1
        if((v(i) >= x(j)) && (v(i) <= x(j+1)))
          s = v(i) - x(j);
          w(i) = y(j) + s * (y(j+1)-y(j))/(x(j+1)-x(j)); 
        endif
      endfor
     endfor
end
