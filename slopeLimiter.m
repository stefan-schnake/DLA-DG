function [f] = slopeLimiter(x,k,f,M,g)
%Applies slope limiter to f.

if k == 0 %Do nothing
    return 
end

%Only works for quadratics or less
assert(k <= 2,'Limiter only supported for k <= 2')

nu = sqrt( (x(2)-x(1))/2 );
dx = x(2)-x(1);

%Get transformation matrix from DG basis to limiter basis
%L = (1/nu)*[1/sqrt(2) 0 3*sqrt(10)/4;0 sqrt(6) 0; 0 0 sqrt(10)/3];
L = (1/nu)*diag([1/sqrt(2),sqrt(6)/2,3*sqrt(10)/4]);
%And inverse (lim to DG)
%L_inv = nu*[sqrt(2) 0 -9/(2*sqrt(2));0 1/sqrt(6) 0;0 0 3/sqrt(10)];
L_inv = nu*diag(1./[1/sqrt(2),sqrt(6)/2,3*sqrt(10)/4]);


N = numel(f)/(k+1); %Number of elements

%Scale down for deg.  %Inverse commutes since matrix is triangular.  
L = L(1:(k+1),1:(k+1));
L_inv = L_inv(1:(k+1),1:(k+1));

%Convert from DG to lim
f_curr = L*f(1:(k+1));
f_next = L*f(k+2:2*(k+1));
for i=2:N-1
    
    %Current index
    idx = (i-1)*(k+1) + 1 : i*(k+1);
    
    %Shift
    f_prev = f_curr;
    f_curr = f_next;
    
    %Calc f_next
    f_next = L*f(idx+(k+1));    
    
    if abs(f_curr(2)) < M*dx^2 %Do not limit
        f_lim = f_curr;
    else
        %fprintf('Limiting on cell %d\n',i);
        f_lim = zeros(size(f_curr));
        f_lim(1) = f_curr(1); %Set constant coeff
        f_lim(2) = minmod(f_curr(2),g*(f_next(1)-f_curr(1)),g*(f_curr(1)-f_prev(1)));  %Linear coeff set by minmod
        %%% f_lim(3) = 0 if quadratic
    end
    
    %Map back to DG space
    f(idx) = L_inv*f_lim;
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

