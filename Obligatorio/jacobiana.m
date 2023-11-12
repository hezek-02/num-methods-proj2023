
function retval = jacobiana(x,S,minCuad)
  c = 299792.458; %km/s
  if (size(x,2) == 3)  
    JF  = [ 
      (x(1)-S(1,1))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),      (x(2)-S(1,2))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),        (x(3)-S(1,3))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2);                    
      (x(1)-S(2,1))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),      (x(2)-S(2,2))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),        (x(3)-S(2,3))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2);
      (x(1)-S(3,1))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),      (x(2)-S(3,2))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),        (x(3)-S(3,3))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2);
      (x(1)-S(4,1))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),      (x(2)-S(4,2))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),        (x(3)-S(4,3))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2);
    ];
  elseif (size(x,2) == 4 && minCuad==0)  
    JF  = [ 
      (x(1)-S(1,1))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),      (x(2)-S(1,2))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),        (x(3)-S(1,3))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),   c ;                        
      (x(1)-S(2,1))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),      (x(2)-S(2,2))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),        (x(3)-S(2,3))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),   c ;
      (x(1)-S(3,1))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),      (x(2)-S(3,2))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),        (x(3)-S(3,3))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),   c ;
      (x(1)-S(4,1))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),      (x(2)-S(4,2))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),        (x(3)-S(4,3))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),   c ;
    ];
   elseif (size(x,2) ==4 && minCuad==8)
    JF  = [
      (x(1)-S(1,1))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),      (x(2)-S(1,2))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2),        (x(3)-S(1,3))/sqrt((x(1) - S(1,1)).^2 + (x(2) - S(1,2)).^2 + (x(3) - S(1,3)).^2) ,   c ;
      (x(1)-S(2,1))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),      (x(2)-S(2,2))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2),        (x(3)-S(2,3))/sqrt((x(1) - S(2,1)).^2 + (x(2) - S(2,2)).^2 + (x(3) - S(2,3)).^2) ,   c ;
      (x(1)-S(3,1))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),      (x(2)-S(3,2))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2),        (x(3)-S(3,3))/sqrt((x(1) - S(3,1)).^2 + (x(2) - S(3,2)).^2 + (x(3) - S(3,3)).^2) ,   c ;
      (x(1)-S(4,1))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),      (x(2)-S(4,2))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2),        (x(3)-S(4,3))/sqrt((x(1) - S(4,1)).^2 + (x(2) - S(4,2)).^2 + (x(3) - S(4,3)).^2) ,   c ;
      (x(1)-S(5,1))/sqrt((x(1) - S(5,1)).^2 + (x(2) - S(5,2)).^2 + (x(3) - S(5,3)).^2),      (x(2)-S(5,2))/sqrt((x(1) - S(5,1)).^2 + (x(2) - S(5,2)).^2 + (x(3) - S(5,3)).^2),        (x(3)-S(5,3))/sqrt((x(1) - S(5,1)).^2 + (x(2) - S(5,2)).^2 + (x(3) - S(5,3)).^2) ,   c ;
      (x(1)-S(6,1))/sqrt((x(1) - S(6,1)).^2 + (x(2) - S(6,2)).^2 + (x(3) - S(6,3)).^2),      (x(2)-S(6,2))/sqrt((x(1) - S(6,1)).^2 + (x(2) - S(6,2)).^2 + (x(3) - S(6,3)).^2),        (x(3)-S(6,3))/sqrt((x(1) - S(6,1)).^2 + (x(2) - S(6,2)).^2 + (x(3) - S(6,3)).^2) ,   c ;
      (x(1)-S(7,1))/sqrt((x(1) - S(7,1)).^2 + (x(2) - S(7,2)).^2 + (x(3) - S(7,3)).^2),      (x(2)-S(7,2))/sqrt((x(1) - S(7,1)).^2 + (x(2) - S(7,2)).^2 + (x(3) - S(7,3)).^2),        (x(3)-S(7,3))/sqrt((x(1) - S(7,1)).^2 + (x(2) - S(7,2)).^2 + (x(3) - S(7,3)).^2) ,   c ;
      (x(1)-S(8,1))/sqrt((x(1) - S(8,1)).^2 + (x(2) - S(8,2)).^2 + (x(3) - S(8,3)).^2),      (x(2)-S(8,2))/sqrt((x(1) - S(8,1)).^2 + (x(2) - S(8,2)).^2 + (x(3) - S(8,3)).^2),        (x(3)-S(8,3))/sqrt((x(1) - S(8,1)).^2 + (x(2) - S(8,2)).^2 + (x(3) - S(8,3)).^2) ,   c ;
    ];

  endif
  retval =JF;
  return;
endfunction
