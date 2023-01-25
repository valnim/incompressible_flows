function x = GaussSeidel(A, b, x0, tol, maxit)
    % A is the coefficient matrix
    % b is the right-hand side vector
    % x0 is the initial guess for the solution
    % tol is the tolerance for the residual
    % maxit is the maximum number of iterations
    n = size(A,1);
    x = x0;
    for k = 1:maxit
        x_prev = x;
        for i = 1:n
            s = 0;
            for j = 1:i-1
                s = s + A(i,j)*x(j);
            end
            for j = i+1:n
                s = s + A(i,j)*x_prev(j);
            end
            x(i) = (b(i) - s)/A(i,i);
        end
        if norm(x-x_prev) < tol
            break;
        end
    end
end