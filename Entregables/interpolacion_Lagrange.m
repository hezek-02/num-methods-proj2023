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
## @deftypefn {} {@var{retval} =} interpolacion_Lagrange (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ezequiel Capó <hezek@heze-Makina>
## Created: 2023-09-22

function p = interpolacion_Lagrange(x, y,v)
    n = length(x) ;%grado del polinomio n-1
    p = zeros(size(v));
    for k = 1:n
        L = ones(size(v));
        for j = 1:n
            if j ~= k
                L = L .*( (v - x(j))./ (x(k) - x(j)));
            end
        end
        p = p + y(k) * L;
    end
end

