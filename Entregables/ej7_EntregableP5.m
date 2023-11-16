
A =@(delta) [1 1 1; delta 0 0; 0 delta 0; 0 0 delta];
y = [1; 0; 0; 0];

delta = [ 1, 0.001,  0.00000001,  eps, sqrt(eps), eps/2, 0.00000000000000000000000001  ];
fprintf('\nSe computará la solución de mínimos cuadrados del sistema Ax=y \n\n')

###Método ec normales
A_cuadrada =@(delta) transpose(A(delta)) * A(delta);
y_fit =@(delta)  transpose(A(delta))*y;

for d=1:length(delta)
  sol_Ec_normal = A_cuadrada(d) \ y_fit(d);
  fprintf('\nSolución con ec. normal con delta = %d, y vector solución: \n',delta(d));
  disp(sol_Ec_normal);
endfor  

###Método HouseHolder
for d=1:3
  A_n = A(d);
  u1 =  A_n(:,1) - norm(A_n(:,1)) * [1;0;0;0];
  rho = 2/(norm(u1,2)^2);
  H1 = eye(4) - rho * u1 * u1';
  A2 = H1*A_n;
  tA2 = A2(2:end,2:end);

  u1 =  tA2(:,1) - norm(tA2(:,1)) * [1;0;0];
  rho = 2/(norm(u1,2)^2);
  tH2 = eye(3) - rho * u1 * u1';
  H2 = eye(4);
  H2(2:end,2:end) = tH2;
  A3 = H2 * A2; 
  tA3 = A3(3:end,3:end);
 
  u1 =  tA3(:,1) - norm(tA3(:,1)) * [1;0];
  rho = 2/(norm(u1,2)^2);
  tH3 = eye(2) - rho * u1 * u1';
  H3 = eye(4);
  H3(3:end,3:end) = tH3;
  A4 = H3 * A3; 
  z = (H1*H2*H3)' * y;
  solHouseHold = A4(1:3, 1:3) \ z(1:3);
  fprintf('\nSolución con HouseHolder, con delta = %d, y vector solución: \n',delta(d));
  disp(solHouseHold);
endfor


###Método SVD
for d=1:length(delta)
  [U, S, V] = svd(A(d));
  z = U'*y;
  w = S(1:3, 1:3) \ z(1:3);
  solSvd = V*w;
  fprintf('\nSolución con SVD, con delta = %d, y vector solución: \n',delta(d));
  disp(solSvd);
endfor  
