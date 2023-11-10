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
## @deftypefn {} {@var{retval} =} gauss_newton (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ezequiel Capó <hezek@heze-Makina>
## Created: 2023-11-10

function [x,k] = gauss_newton(x0,y, JF, F, tol, itMax)
  x = x0;
  k = 0;
  dif = inf;
  while (k<itMax) && (dif > tol)
      r = y' - F(x);
      [U,S,V] = svd(JF(x));
      z = U'*((-1)*r);
      w = S(1:4,1:4) \ z;
      dif = V*w;
      x = x + dif;
      k = k+1;
  endwhile
endfunction
