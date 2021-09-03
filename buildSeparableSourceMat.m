function [SC] = buildSeparableSourceMat(x,v,k,t,fcell)
%Builds (f(x,v),phi_j)_{\W} where f(x,v)=fx(x)fv(v)

%Here fcell is a cell with 3 components
%fcell{i,1} = function of x,t
%fcell{i,2} = function of t
%fcell{i,3} = funciton of v,t

SC = cell(size(fcell));

num_x = numel(x)-1;
num_v = numel(v)-1;

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

NN = size(fcell,1);
for l=1:NN
    SC{l,1} = zeros((k+1)*num_x,1);
    SC{l,3} = zeros((k+1)*num_v,1);
end

%%%Build x vectors
[leg_vals,~,~,~] = buildLegendre(10,k);
leg_vals = leg_vals/sqrt(jac_x);
test_ref = repmat(w_ref',k+1,1).*leg_vals; %Weights included with the test functions

for i=1:num_x
    quad_x = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
    for l=1:NN
        SC{l,1}((k+1)*(i-1)+1:(k+1)*i) = fcell{l,1}(quad_x,t)*test_ref'*jac_x;
    end
end


%Build v vectors
[leg_vals,~,~,~] = buildLegendre(10,k);
leg_vals = leg_vals/sqrt(jac_v);
test_ref = repmat(w_ref',k+1,1).*leg_vals; %Weights included with the test functions

for i=1:num_v
    quad_v = quad_ref*(v(i+1)-v(i))/2 + (v(i+1)+v(i))/2;
    for l=1:NN
        SC{l,3}((k+1)*(i-1)+1:(k+1)*i) = fcell{l,3}(quad_v,t)*test_ref'*jac_v;
    end
end

%Normalize and add time scalar
for l=1:NN
    norm_x = norm(SC{l,1},2);
    SC{l,1} = SC{l,1}/norm_x;
    norm_v = norm(SC{l,3},2);
    SC{l,3} = SC{l,3}/norm_v;
    SC{l,2} = fcell{l,2}(t)*norm_x*norm_v;
end



end

