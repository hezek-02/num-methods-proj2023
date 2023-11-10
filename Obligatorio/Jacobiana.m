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
## @deftypefn {} {@var{retval} =} Jacobiana(@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ezequiel Capó <hezek@heze-Makina>
## Created: 2023-11-09

function retval = Jacobiana(x,S)
  if (size(x,2) == 3)  
    Jf = [ 
      (x(1)-S(1,1))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),      (x(2)-S(1,2))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),        (x(3)-S(1,3))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2);                    
      (x(1)-S(2,1))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),      (x(2)-S(2,2))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),        (x(3)-S(2,3))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2);
      (x(1)-S(3,1))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),      (x(2)-S(3,2))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),        (x(3)-S(3,3))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2);
      (x(1)-S(4,1))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),      (x(2)-S(4,2))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),        (x(3)-S(4,3))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2);
    ];
  elseif (size(x,2) == 4)  
    c = 299792.458; %km/s
    Jf = [ 
      (x(1)-S(1,1))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),      (x(2)-S(1,2))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),        (x(3)-S(1,3))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),   c ;                        
      (x(1)-S(2,1))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),      (x(2)-S(2,2))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),        (x(3)-S(2,3))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),   c ;
      (x(1)-S(3,1))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),      (x(2)-S(3,2))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),        (x(3)-S(3,3))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),   c ;
      (x(1)-S(4,1))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),      (x(2)-S(4,2))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),        (x(3)-S(4,3))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),   c ;
    ];
  end
  retval = Jf;
  return;
endfunction