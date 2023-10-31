
function retval = newton_raphson (x0, f, Jf, tol, max_it)
    x = x0;  % Estimación inicial
    k = 0;  % Inicializar el contador de iteraciones

    while k < max_it
        b = -f(x);
        A = Jf(x);

        % Resolver el sistema lineal Ad = b para encontrar d
        d = A \ b;

        x = x + d;  % Actualizar x

        % Comprobar el criterio de parada
        if norm(d) < tol
            break;  % Detenerse si se cumple el criterio de parada
        end

        k = k + 1;  % Incrementar el contador de iteraciones
    end


endfunction
