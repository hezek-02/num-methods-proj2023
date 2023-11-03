
function [x,k] = newton_raphson (x0, Jf, F, itMax, critExtra)
  k = 0;
  x = x0;
  xAnt = x +100.* 100;
  if critExtra == 1
    while norm(x-xAnt,2)  > 1e-6 &&  k < itMax
      xAnt = x;
      A = Jf(x);
      d = A\-F(x);
      x = x + d'; %tomo la trasp, x es 1xn y d es nx1
      k = k+1;
    endwhile
  else
      while  k < itMax
      xAnt = x;
      A = Jf(x);
      d = A\-F(x);
      x = x + d'; %tomo la trasp, x es 1xn y d es nx1
      k = k+1;
    endwhile
  endif
  return;
endfunction
