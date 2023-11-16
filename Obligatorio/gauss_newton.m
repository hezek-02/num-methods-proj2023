function [x,k] = gauss_newton(x0, y, JF, F, tol, itMax)
  x = x0';
  k = 0;
  dif = inf;
  while (k <itMax) && (norm(dif) > tol)
      r = F(x') - y;
      [U,S,V] = svd(JF(x'));
      z = U' * (-r);
      n = size(x,1);
      w = S(1:n,1:n) \ z(1:n);
      dif = V*w;
      x = x + dif;
      k = k+1;
 endwhile
endfunction
