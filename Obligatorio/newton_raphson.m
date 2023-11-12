
function [x,k] = newton_raphson (x0, Jf, F, tol, itMax)
  k = 0;
  x = x0;
  xAnt = inf;
    while norm(x-xAnt,2)  > tol &&  k < itMax
      xAnt = x;
      A = Jf(x);
      d = A\-F(x);
      x = x + d'; %tomo la trasp, x es 1xn y d es nx1
      k = k+1;
    endwhile
  return;
endfunction
