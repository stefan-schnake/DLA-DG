function [Lf] = wenoLimiter(x,k,f,M)
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

Lf = f;


N = numel(f)/(k+1); %Number of elements

%Scale down for deg.  %Inverse commutes since matrix is triangular.  
L = L(1:(k+1),1:(k+1));
L_inv = L_inv(1:(k+1),1:(k+1));

%WENO Parameters
gamma = [0.001,0.998,0.001];
epsilon = 1e-6;
r = 2;

%Convert from DG to lim
for i=2:N-1
    
    %Current index
    idx = (i-1)*(k+1) + 1 : i*(k+1);
    
    %Shift
    left = L*f(idx-(k+1));
    curr = L*f(idx);
    rght = L*f(idx+(k+1));
       
    
    if abs(curr(2)) < M*dx^2 %Do not limit
        alpha = curr(2);
    else
        alpha = minmod(curr(2),rght(1)-curr(1),curr(1)-left(1));  %Linear coeff set by minmod
    end
    if abs(alpha-curr(2)) > 1e-14 %%Need to build WENO parameters
        if k == 1
            beta = 4*[rght(2)^2,curr(2)^2,left(2)^2];
        else
            beta = [4*rght(2)^2+32*rght(2)*rght(3)+(16*13/3+64)*rght(3)^2,...
                    4*curr(2)^2+(16/3+64)*curr(3)^2,...
                    4*left(2)^2-32*left(2)*left(3)+(16*13/3+64)*left(3)^2];
        end
        omega_bar = gamma./(beta+epsilon).^r;
        omega = omega_bar/sum(omega_bar(:));
        
        %    Need to normal p_0 (rght) on I_{i-1} and
        % p_2 (left) on I_{i+1} so that their integrals of I_i are the same
        % as the average value of p_1 on I_i;  
        if k == 1
            rght(1) = rght(1) - dx*[1  2]*rght + curr(1);
            left(1) = left(1) - dx*[1 -2]*left + curr(1);
        else
            rght(1) = rght(1) - dx*[1  2 14/3]*rght + curr(1);
            left(1) = left(1) - dx*[1 -2 14/3]*left + curr(1);
        end
        
        %Redefine curr
        curr = omega(1)*rght + omega(2)*curr + omega(3)*left;
    end
    
    %Map back to DG space
    Lf(idx) = L_inv*curr;
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