function [ru] = slopeLimiter2D(x,v,k,u,M)
%2D slope limiter for linear DG functions

if k == 0
    ru = u;
    return 
end

assert(k < 2,'Error: Slope limiter only works for linear and constant polynomials')

quiet = true;

ru = u;

num_x = numel(x)-1;
num_v = numel(v)-1;

dx = (x(2)-x(1));
dv = (v(2)-v(1));

%Transformation from DG basis to Limiter basis
L = diag([sqrt(2/(dx*dv)),6/sqrt(dx*dv),6/sqrt(dx*dv),6/(2*sqrt(dx*dv))]);
%Transformation from Limiter basis to DG basis
Linv = inv(L);

for i=2:num_x-1
    for j=2:num_v-1
        %Convert to proper basis
        %Get local index
        start = ((i-1)*num_v+(j-1))*(k+1)^2;
        curr = L*u(start+1:start+(k+1)^2);
        
        %Get elements to up, down, left, and right of u
        idx = start + (k+1)^2;
        up = L*u(idx+1:idx+(k+1)^2);
        
        idx = start - (k+1)^2;
        dn = L*u(idx+1:idx+(k+1)^2);
        
        idx = start - num_v*(k+1)^2;
        lt = L*u(idx+1:idx+(k+1)^2);
        
        idx = start + num_v*(k+1)^2;
        rt = L*u(idx+1:idx+(k+1)^2);
        
        if abs(curr(2)) >= M*dv^2
            curr(2) = minmod(curr(2),up(1)-curr(1),curr(1)-dn(1));
            curr(4) = 0;
            if ~quiet
                fprintf('Limit in v for cell %d\n',(i-1)*num_v+j);
            end
        end
        if abs(curr(3)) >= M*dx^2
            curr(3) = minmod(curr(3),rt(1)-curr(1),curr(1)-lt(1));
            curr(4) = 0;
            if ~quiet
                fprintf('Limit in x for cell %d\n',(i-1)*num_v+j);
            end
        end
        
        %Map back and store
        ru(start+1:start+(k+1)^2) = Linv*curr;
    end
end


end

function z = minmod(a,b,c)
    if abs(sign(a)+sign(b)+sign(c)) == 3 %Equiv to sign(a) = sign(b) = sign(c)
        s = sign(a);
        z = s*min(abs([a,b,c]));
    else
        z = 0;
    end
end

