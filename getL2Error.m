function [z] = getL2Error(x,k,u,f)
%Outputs \|u-f\|_{L^2(\W)}

num_x = numel(x)-1;
jac_x = (x(2)-x(1))/2;
z = 0;

[quad_x,w] = lgwt(12,-1,1);
[leg_vals,~,~,~] = buildLegendre(12,k);
leg_vals = leg_vals'/sqrt(jac_x);
for i=1:num_x
    aff_x = (x(i)*(1-quad_x)+x(i+1)*(1+quad_x))/2;   
    F = f(aff_x);
    U = leg_vals*u((i-1)*(k+1)+1:i*(k+1));
    Z = (F-U);
    z = z + (Z'*(w.*Z))*jac_x;
end

z = sqrt(z);
end

