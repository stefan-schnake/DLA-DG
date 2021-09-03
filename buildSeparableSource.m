function [F] = buildSeparableSource(x,v,k,fx,fv)
%Builds (f(x,v),phi_j)_{\W} where f(x,v)=fx(x)fv(v)

num_x = numel(x)-1;
num_v = numel(v)-1;

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

[leg_vals,~,~,~] = buildLegendre(10,k);
leg_vals = leg_vals/sqrt(jac_x);
test_ref = repmat(w_ref',k+1,1).*leg_vals; %Weights included with the test functions

block_x = zeros(k+1,num_x);
for i=1:num_x
    quad_x = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
    block_x(:,i) = fx(quad_x)*test_ref'*jac_x;
end

[leg_vals,~,~,~] = buildLegendre(10,k);
leg_vals = leg_vals/sqrt(jac_v);
test_ref = repmat(w_ref',k+1,1).*leg_vals; %Weights included with the test functions

block_v = zeros(k+1,num_v);
for i=1:num_v
    quad_v = quad_ref*(v(i+1)-v(i))/2 + (v(i+1)+v(i))/2;
    block_v(:,i) = fv(quad_v)*test_ref'*jac_v;
end

F = zeros((k+1)^2*num_x*num_v,1);
count = 1;
for i=1:num_x
    for j=1:num_v
        F(count:count+(k+1)^2-1) = kron(block_x(:,i),block_v(:,j));
        count = count+(k+1)^2;
    end
end

end

