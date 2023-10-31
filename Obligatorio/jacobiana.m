## Copyright (C) 2023 facundo
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
## @deftypefn {} {@var{retval} =} jacobiana (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: facundo <ususario@H510M-H>
## Created: 2023-10-31

function retval = jacobiana (f, x0, epsilon)
    n = length(f(x0));
    m = length(x0);
    J = zeros(n, m);

    for i = 1:n
        for j = 1:m
            x_plus_epsilon = x0;
            x_plus_epsilon(j) = x0(j) + epsilon;
            J(i, j) = (f(x_plus_epsilon)(i) - f(x0)(i)) / epsilon;
        end
    end

endfunction
