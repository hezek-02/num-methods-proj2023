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
## @deftypefn {} {@var{retval} =} perturbaciones (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ezequiel Capó <hezek@heze-Makina>
## Created: 2023-11-10

function retval = perturbaciones (tipo)
if (tipo == 1)
    pert = zeros(16,4);
      for j = 0:15
        tira = dec2bin(j,4);
        for k = 1:4 
            pert(j+1,k) = str2num(tira(k));
        end
    end
  pert = 2 * -10e-8 * pert + 10e-8; 
 elseif (tipo == 2)
    pert = zeros(256,8);
    for j = 0:255
        tira = dec2bin(j,8);
        for k = 1:8 
            pert(j+1,k) = str2num(tira(k));
        end
    end
    pert = 2*pert-1;
endif
retval = pert;
return;
endfunction

%Para tipo=1
                  %t1       %t2         %t3          %t4
##delta_ti = [ 10e-8,    10e-8,    10e-8,     10e-8;
##                 -10e-8,   10e-8,    10e-8,      10e-8;
##                  10e-8,    -10e-8,   10e-8,      10e-8;
##                  10e-8,    10e-8,    -10e-8,     10e-8;
##                  10e-8,    10e-8,     10e-8,    -10e-8;
##                  -10e-8,   -10e-8,   10e-8,      10e-8;
##                  -10e-8,   10e-8,    -10e-8,     10e-8;
##                  -10e-8,   10e-8,    10e-8,     -10e-8;
##                  10e-8,    -10e-8,   -10e-8,     10e-8;
##                  10e-8,    -10e-8,    10e-8,    -10e-8;
##                  10e-8,    10e-8,    -10e-8,    -10e-8;
##                  -10e-8,   -10e-8,   -10e-8,    10e-8;
##                  -10e-8,   10e-8,    -10e-8,    -10e-8;
##                  -10e-8,   -10e-8,   10e-8,     -10e-8;
##                  10e-8,    -10e-8,   -10e-8,    -10e-8;
##                  -10e-8,   -10e-8,   -10e-8,    -10e-8];
