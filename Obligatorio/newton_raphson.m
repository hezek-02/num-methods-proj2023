
function x = newton_raphson (x0, Jf, F, itMax)
  x = x0;
  k = 0;
  while k < itMax
    A = Jf(x(1),x(2),x(3),x(4));
    d = A\-F(x(1),x(2),x(3),x(4));
    x = x + d'; %tomo la trasp, x es 1x4 y d 4x1
    k = k+1;
  endwhile
  return;
endfunction
