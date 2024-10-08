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
## @deftypefn {} {@var{retval} =} aproximarRaiz (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ezequiel Capó <hezek@heze-Makina>
## Created: 2023-08-01

function aprRaiz = aproximarRaiz (it,raizCalcular)
  #raiz entre raizCalcular/ao y ao
  if  it == 0
    disp("Error");
    aprRaiz = 0;
    return
  else
    while (abs(it*it - raizCalcular) > 1e-6)
      it = 0.5*(it + raizCalcular/it);
    endwhile
  endif
  aprRaiz = it;
endfunction
















