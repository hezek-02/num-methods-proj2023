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
## @deftypefn {} {@var{retval} =} combinaciones_delta_ti (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ezequiel Capó <hezek@heze-Makina>
## Created: 2023-11-09

function combinaciones = combinaciones_delta_ti()
    valores = [1e-8, 1e8];
    num_valores = 2;
    num_espacios = 4;
    % Calcula el número total de combinaciones posibles
    num_combinaciones = num_espacios^num_valores;
    % Inicializa la matriz para almacenar las combinaciones
    combinaciones = zeros(num_combinaciones, num_espacios);

    % Genera todas las combinaciones
    for i = 1:num_combinaciones
        indice_actual = i;
        for j = 1:num_espacios
            combinaciones(i, j) = valores(mod(indice_actual-1, num_valores) + 1);
            indice_actual = floor((indice_actual-1) / num_valores);
        end
    end

endfunction
