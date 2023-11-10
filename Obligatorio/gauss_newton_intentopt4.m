




function [x,k] = gauss_newton(x0,y,JF,F,tol,itMax)
  x = x0;
  k = 0;
  dif = inf;
  while (k<itMax) && (dif > tol)
      r = y' - F
      [U,S,V] = svd(JF);
      z = U'*((-1)*r);
      w = S(1:3,1:3) \ z;
      dif = V*w;
      x = x + dif
      k = k+1;
  endwhile
