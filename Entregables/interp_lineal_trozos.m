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
## @deftypefn {} {@var{retval} =} interp_lineal_trozos (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ezequiel Capó <hezek@heze-Makina>
## Created: 2023-10-02

function w = interp_lineal_trozos(x, y, v)
    delta_y = diff(y);
    delta_x = diff(x);

    delta_dividido = delta_y ./ delta_x;
    delta_dividido(length(x))=delta_dividido(length(delta_dividido));
    % Encontrar el índice i tal que x(i) <= v < x(i+1)
    i = 1;
    while x(i) < v 
      i = i + 1;
    endwhile
    disp(i);
    % Calcular la variable local s y evaluar L(v)
    s = v - x(i);
    w = y(i) + s * delta_dividido(i);
end