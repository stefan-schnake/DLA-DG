function [F] = buildSeparableVOnlySource(x,v,k,fv,g)
%Builds F_j(x) := (fv(v),g(x,v))_{\W_v}

num_x = numel(x)-1;
num_v = numel(v)-1;

jac_v = (v(2)-v(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

[leg_vals,~,~,~] = buildLegendre(10,k);
leg_vals = leg_vals/sqrt(jac_v);
test_ref = repmat(w_ref',k+1,1).*leg_vals; %Weights included with the test functions


block_v = zeros(k+1,num_v);
for i=1:num_v
    quad_v = quad_ref*(v(i+1)-v(i))/2 + (v(i+1)+v(i))/2;
    block_v(:,i) = fv(quad_v)*test_ref'*jac_v;
end

%F = zeros((k+1)*num_x,(k+1)^2*num_x*num_v);
F = zeros((k+1)*num_x,1);
count = 1;
for i=1:num_x
    for j=1:num_v
        for l=1:k+1
            F((i-1)*(k+1)+l) = F((i-1)*(k+1)+l) + ...
                    block_v(:,j)'*g(count:count+(k+1)-1);
            count = count + (k+1);
        end
    end
end

end

