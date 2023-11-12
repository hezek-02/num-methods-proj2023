
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
  pert = 2 * -10e-8 * pert + 10e-8; 
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
