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
vector = [1 2 3 4 5 6];
normaVec = norma(vector);
matriz = [1 2 3 4 5 6; 3 2 1 3 2 1];
#disp("norma euclidea personal: " + normaVec);
#disp("norma euclidea octave: " + norm(vector));
disp(normaVec);
disp(norm(vector));
raizProx = aproximarRaiz(6,25);
disp(raizProx);
disp(sustitucion(matriz,vector,3));
toc