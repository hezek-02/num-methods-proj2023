
%Entregable ej11-Practico2 Matrices de Hilbert

%a
n1=5;
n2=10;

xn1 = ones(n1,1); %ya esta traspuesta
xn2 = ones(n2,1);

Hn1 = hilb(n1);
Hn2 = hilb(n2);

bn1 = Hn1 * xn1;
bn2 = Hn2 * xn2;

sol_n1 = Hn1 \ bn1;
sol_n2 = Hn2 \ bn2;

disp("Solución \ Hilbert 5x5: ");
disp(sol_n1);
disp("\n Solución \ Hilbert 10x10: ");
disp(sol_n2);
%b
%errores
disp("ERRORES:");
disp(" ");

en1 = xn1 - sol_n1;
en2 = xn2 - sol_n2;

norma_en1_1 = norm(en1,1);
norma_en1_2 = norm(en1,2);
norma_en1_inf = norm(en1,inf);
disp("norma del error al operar con Hilbert 5x5");
disp("norma lp=1:");
disp(norma_en1_1);
disp(" ");
disp("norma lp=2:");
disp(norma_en1_2);
disp(" ");
disp("norma Infinito:");
disp(norma_en1_inf);
disp(" ");

norma_en2_1 = norm(en2,1);
norma_en2_2 = norm(en2,2);
norma_en2_inf = norm(en2,inf);
disp("norma del error al operar con Hilbert 10x10");
disp("norma lp=1:");
disp(norma_en2_1);
disp(" ");
disp("norma lp=2:");
disp(norma_en2_2);
disp(" ");
disp("norma Infinito:");
disp(norma_en2_inf);
disp(" ");

%residuos
disp("RESIDUOS:");
disp(" ");
rn1 = bn1 - Hn1*sol_n1;
rn2 = bn2 - Hn2*sol_n2;

norma_rn1_1 = norm(rn1,1);
norma_rn1_2 = norm(rn1,2);
norma_rn1_inf = norm(rn1,inf);
disp("norma del residuo al operar con Hilbert 5x5");
disp("norma lp=1:");
disp(norma_rn1_1);
disp(" ");
disp("norma lp=2:");
disp(norma_rn1_2);
disp(" ");
disp("norma Infinito:");
disp(norma_rn1_inf);
disp(" ");

norma_rn2_1 = norm(rn2,1);
norma_rn2_2 = norm(rn2,2);
norma_rn2_inf = norm(rn2,inf);
disp("norma del residuo al operar con Hilbert 10x10");
disp("norma lp=1:");
disp(norma_rn2_1);
disp(" ");
disp("norma lp=2:");
disp(norma_rn2_2);
disp(" ");
disp("norma Infinito:");
disp(norma_rn2_inf);
disp(" ");

%c
disp("NÚMEROS DE CONDICIÓN:");
disp(" ");
disp("Número de condicion al operar con Hilbert 5x5 (con cond):");
disp("con norma lp=1");
disp(cond(Hn1,1));
disp(" ");
disp("con norma lp=2");
disp(cond(Hn1,2));
disp(" ");
disp("con norma infinito");
disp(cond(Hn1,inf));
disp(" ");

disp("Número de condicion al operar con Hilbert 10x10 (con cond):");
disp("con norma lp=1");
disp(cond(Hn2,1));
disp(" ");
disp("con norma lp=2");
disp(cond(Hn2,2));
disp(" ");
disp("con norma infinito");
disp(cond(Hn2,inf));
disp(" ");


disp(" ");
disp("Número de condicion al operar con Hilbert 5x5:");
disp("con norma lp=1");
disp(norm(Hn1,1)*norm(invhilb(n1),1));
disp(" ");
disp("con norma lp=2");
disp(norm(Hn1,2)*norm(invhilb(n1),2));
disp(" ");
disp("con norma infinito");
disp(norm(Hn1,inf)*norm(invhilb(n1),inf));
disp(" ");

disp("Número de condicion al operar con Hilbert 10x10:");
disp("con norma lp=1");
disp(norm(Hn2,1)*norm(invhilb(n2),1));
disp(" ");
disp("con norma lp=2");
disp(norm(Hn2,2)*norm(invhilb(n2),2));
disp(" ");
disp("con norma infinito");
disp(norm(Hn2,inf)*norm(invhilb(n2),inf));
disp(" ");






