

function [x,k] = gauss_newton(x0, y, JF, F, tol, itMax)
  x = x0;
  k = 0;
  dif = inf;
  while (k <itMax) && (norm(dif) > tol)
      r = F(x) - y;
      [U,S,V] = svd(JF(x));
      z = U' * (-r);
      w = S(1:4,1:4) \ z(1:4);
      dif = V*w;
      x = x + dif;
      k = k+1;
 endwhile
endfunction
