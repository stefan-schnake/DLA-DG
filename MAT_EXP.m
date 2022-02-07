function z = MAT_EXP(A,t,z,tol)
    %Matrix exponential
    r = z;
    for i=1:20
        r = -(A*r);
        z = z + 1/factorial(i)*t^i*r;
        if 1/factorial(i)*t^i*norm(r) < tol
            break
        end
    end
end