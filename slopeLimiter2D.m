function [f] = slopeLimiter2D(x,v,k,f,M,g)
%Applies slope limiter to f.

if nargin < 6
    g = 1;
end

quiet = false;

%Only works for quadratics or less
assert(k <= 2,'Limiter only supported for k <= 2')

num_x = numel(x)-1;
num_v = numel(v)-1;

if k == 0 %Do nothing
    return 
elseif k == 1
    lim_deg = [1,2,3];
    zero_deg = [4];
else
    lim_deg = [1,2,4];
    zero_deg = [3,5,6,7,8,9];
end

nu_x = sqrt( (x(2)-x(1))/2 );
nu_v = sqrt( (v(2)-v(1))/2 );
dx = x(2)-x(1);
dv = v(2)-v(1);

%Get transformation matrix from DG basis to limiter basis
%L = (1/nu)*[1/sqrt(2) 0 3*sqrt(10)/4;0 sqrt(6) 0; 0 0 sqrt(10)/3];
L = diag([1/sqrt(2),sqrt(6)/2,3*sqrt(10)/4]);
%And inverse (lim to DG)
%L_inv = [sqrt(2) 0 -9/(2*sqrt(2));0 1/sqrt(6) 0;0 0 3/sqrt(10)];
L_inv = diag(1./[1/sqrt(2),sqrt(6)/2,3*sqrt(10)/4]);

%Scale down for deg and kron.  %Inverse commutes since matrix is triangular.  
L = kron(L(1:(k+1),1:(k+1)),L(1:(k+1),1:(k+1)))/(nu_x*nu_v); 
    L = L(lim_deg,lim_deg);
L_inv = kron(L_inv(1:(k+1),1:(k+1)),L_inv(1:(k+1),1:(k+1)))*(nu_x*nu_v);
    L_inv = L_inv(lim_deg,lim_deg);

%Convert from DG to lim
for i=1:num_x
    for j=1:num_v
    
        %Current index
        ele_idx = (i-1)*num_v + j;
        idx = (ele_idx-1)*(k+1)^2 + 1 : ele_idx*(k+1)^2;
        
        f_curr = L*f(idx(lim_deg));
        zero_out = false; %Determines if we're zeroing out the 
                            %no needed coefficients
        
        if (i > 1) && (i < num_x) %Limit in x
            if abs(f_curr(3)) >= M*dx^2
                %Get avg values to left and right
                avg_left = f( ((i-2)*num_v+j)*(k+1)^2 );
                avg_rght = f( ( i   *num_v+j)*(k+1)^2 );
                f_curr(3) = minmod(f_curr(3),g*(avg_left-f_curr(1)),g*(f_curr(1)-avg_rght));
                zero_out = true;
            end
        end

        if (j > 1) && (j < num_v) %Limit in x
            if abs(f_curr(2)) >= M*dv^2
                %Get avg values to top and bottom
                avg_up = f( ((i-1)*num_v+j-1)*(k+1)^2 );
                avg_dn = f( ((i-1)*num_v+j+1)*(k+1)^2 );
                f_curr(2) = minmod(f_curr(2),g*(avg_up-f_curr(1)),g*(f_curr(1)-avg_dn));
                zero_out = true;
            end
        end
        
        if zero_out %We've limited
            if quiet
                fprintf('Limiting element %d\n',ele_idx);
            end
            f(idx(lim_deg)) = L_inv*f_curr;
            f(idx(zero_deg)) = zeros(size(zero_deg));
        end
    
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

