function [F] = buildSeparableSourceX(x,k,f)
%Builds (f(x),phi_j)_{\W_x}

num_x = numel(x)-1;

jac_x = (x(2)-x(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

[leg_vals,~,~,~] = buildLegendre(10,k);
leg_vals = leg_vals/sqrt(jac_x);
test_ref = repmat(w_ref',k+1,1).*leg_vals; %Weights included with the test functions

block_x = zeros(k+1,num_x);
for i=1:num_x
    quad_x = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
    block_x(:,i) = f(quad_x)*test_ref'*jac_x;
end

F = zeros((k+1)*num_x,1);
count = 1;
for i=1:num_x
    F(count:count+(k+1)-1) = block_x(:,i);
    count = count+(k+1);
end

end

