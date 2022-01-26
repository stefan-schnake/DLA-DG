function [z] = getL2Error(x,v,k,u,utrue)
% Returns ||u_h-u_true||_{L^2)

num_x = numel(x)-1;
num_v = numel(v)-1;
poi = 8;

label_x = 'x';
label_y = 'v'; 

[quad_ref, w_ref]  = lgwt(poi,-1,1);
quad_ref = quad_ref';
%quad_ref = fliplr(quad_ref);

[leg_vals,~,~,~] = buildLegendre(poi,k);
%leg_vals = fliplr(leg_vals);

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

X = zeros(poi*num_x,1);
V = zeros(poi*num_v,1);
for i=1:num_x
    X((i-1)*poi+1:i*poi) = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
end
for i=1:num_v
    V((i-1)*poi+1:i*poi) = quad_ref*(v(i+1)-v(i))/2 + (v(i+1)+v(i))/2;
end

z = 0;
count = 1;
for i=1:num_x
    for j=1:num_v
        coeff = u(count:count+(k+1)^2-1);
        U_h = zeros(poi);
        for ii=1:(k+1)
            for jj=1:(k+1)
                U_h = U_h + coeff((ii-1)*(k+1)+jj)*((leg_vals(jj,:)/sqrt(jac_v))'*(leg_vals(ii,:)/sqrt(jac_x)));
            end
        end
        [XX,VV] = meshgrid(X((i-1)*poi+1:i*poi),V((j-1)*poi+1:j*poi));
        U_true = utrue(XX,VV);
        Z = (U_h-U_true).^2;
        z = z + (w_ref'*Z*w_ref)*jac_x*jac_v;
        count = count + (k+1)^2;
    end
end
z = sqrt(z);

end


