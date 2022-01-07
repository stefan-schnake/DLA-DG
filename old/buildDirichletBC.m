function [F] = buildDirichletBC(x,v,k,u)
%Build Dirichlet BC for (div_x(au),v) where a=(-y,x)

num_x = numel(x)-1;
num_v = numel(v)-1;

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';



F = zeros((k+1)^2*num_x*num_v,1);

[leg_vals,~,leg_edge_vals,~] = buildLegendre(10,k);

%Traverse v edges first
for i=1:num_x 
    quad_x = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
    
    if (x(i+1)+x(i)) < 0 %Inflow on top of domain
        blockstart = (k+1)^2*(i*num_v-1);        
        u_val = u(quad_x,v(end));
        sel = 2;
        nu = 1;
    else %Inflow on bottom of domain
        blockstart = (k+1)^2*((i-1)*num_v);
        u_val = u(quad_x,v(1));
        sel = 1;
        nu = -1;
    end
    
    indices = blockstart+1:blockstart+(k+1)^2;
    
    %Contruct x integral
    temp_integral = (u_val.*w_ref'.*quad_x)*(leg_vals/sqrt(jac_x))'*jac_x;
    
    %Kron with edges in v
    F(indices) = F(indices) + ...
        kron(temp_integral',nu*leg_edge_vals(:,sel)/sqrt(jac_v));
end

for j=1:num_x %Next x edges
    
    quad_v = quad_ref*(v(j+1)-v(j))/2 + (v(j+1)+v(j))/2;
    
    if (v(j+1)+v(j)) < 0 %Inflow on left edge
        blockstart = (k+1)^2*(j-1);
        u_val = u(x(1),quad_v);
        sel = 1;
        nu = -1;
    else %Inflow on right edge
        blockstart = (k+1)^2*((num_x-1)*num_v+j-1);
        u_val = u(x(end),quad_v);
        sel = 2;
        nu = 1;
    end
    
    indices = blockstart+1:blockstart+(k+1)^2;
    
    temp_integral = (u_val.*w_ref'.*(-quad_v))*(leg_vals/sqrt(jac_v))'*jac_v;
    
    F(indices) = F(indices) + ...
        kron(nu*leg_edge_vals(:,sel)/sqrt(jac_x),temp_integral');
end


end