function [BC] = buildDirichletMatBC(x,v,k,u)
%Build Dirichlet BC for (div_x(au),v) where a=(-y,x)

num_x = numel(x)-1;
num_v = numel(v)-1;

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

BC = cell(4,3);
for i=1:4
    BC{i,1} = zeros((k+1)*num_x,1);
    BC{i,2} = 1;
    BC{i,3} = zeros((k+1)*num_v,1);
end

[leg_vals,~,leg_edge_vals,~] = buildLegendre(10,k);

%Traverse v edges first
for i=1:num_x 
    quad_x = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
    
    if (x(i+1)+x(i)) < 0 %Inflow on top of domain
        %non-zero portion is ENE strip of matrix       
        u_val = u(quad_x,v(end));
    else %Inflow on bottom of domain
        %non-zero portion is WSW strip of matrix        
        u_val = u(quad_x,v(1));
    end
    
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Contruct x integral
    temp_integral = (u_val.*w_ref'.*quad_x)*(leg_vals/sqrt(jac_x))'*jac_x;
    
    %Kron with edges in v
    if (x(i+1)+x(i)) < 0
        BC{1,1}(indices) = temp_integral';
        %kron(temp_integral',nu*leg_edge_vals(:,sel)/sqrt(jac_v));
    else
        BC{2,1}(indices) = temp_integral';
    end
end

sel = 2;
nu = 1;
BC{1,3}(end-k:end) = nu*leg_edge_vals(:,sel)'/sqrt(jac_v);
sel = 1;
nu = -1;
BC{2,3}(1:k+1) = nu*leg_edge_vals(:,sel)'/sqrt(jac_v);

%Add v portion
for j=1:num_x %Next x edges
    
    quad_v = quad_ref*(v(j+1)-v(j))/2 + (v(j+1)+v(j))/2;
    
    if (v(j+1)+v(j)) < 0 %Inflow on left edge
        u_val = u(x(1),quad_v);
    else %Inflow on right edge
        u_val = u(x(end),quad_v);
    end
    
    blockstart = (k+1)*(j-1);
    indices = blockstart+1:blockstart+(k+1);
    
    temp_integral = (u_val.*w_ref'.*(-quad_v))*(leg_vals/sqrt(jac_v))'*jac_v;
    if (v(j+1)+v(j)) < 0
        BC{3,3}(indices) = temp_integral';
    else
        BC{4,3}(indices) = temp_integral';
    end
end

sel = 1;
nu = -1;
BC{3,1}(1:k+1) = nu*leg_edge_vals(:,sel)'/sqrt(jac_x);
sel = 2;
nu = 1;    
BC{4,1}(end-k:end) = nu*leg_edge_vals(:,sel)'/sqrt(jac_x);



end
