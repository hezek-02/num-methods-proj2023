
At = [1,0,0,0;2,3,0,0];
A = At';

yt=[1,2,3,4];

y=yt';

xsol=(At*A)\(At*y)
b=At*A
At*y
res = norm(A*xsol-y,2)