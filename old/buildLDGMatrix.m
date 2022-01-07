function [L] = buildLDGMatrix(x,v,k)
%Stiffness matrix for Laplacian -Delta with homogeneous Neumann BCs.
% (\grad g , tau)_\W - <[g],{tau}>_{E_h^I}. 
% g\in V_h, tau=(tau_1,tau_2)\in [V_h]^2

%Each element T\in Th is T=T_x \times T_v.

num_x = numel(x)-1;
num_v = numel(v)-1;

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

[leg_vals,leg_der_vals,leg_edge_vals,~] = buildLegendre(10,k);

block_x = cell(10,num_x);
for i=1:num_x   
    %Create quadrature points on x cell
    quad_x = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
    
    %%%Create appropriate integrals in the x direction.  
    %Each term in this decomposition is seperable so we can compute the
    %integrals in the x-direction and apply them to each v cell
    
    %(phi_i,phi_j)_{T_x}
    block_x{1,i} = (leg_vals/sqrt(jac_x))*(w_ref.*leg_vals'/sqrt(jac_x))*jac_x;
    %(\grad phi_i,phi_j)_{T_x}
    block_x{2,i} = (leg_vals/sqrt(jac_x))*(w_ref.*leg_der_vals'/(jac_x*sqrt(jac_x)))*jac_x;
    
    %Create edge integrals 
    %Note: Since normals are coordinate aligned, nu\cdot\tau kills tau_2.
    block_x{3,i} = zeros(k+1);
    block_x{4,i} = zeros(k+1);
    block_x{5,i} = zeros(k+1);
    if i > 1 %lower edge integral
        %First edge integrals for functions who share same volume
        jump = -leg_edge_vals(:,1)/sqrt(jac_x);
        avg =   leg_edge_vals(:,1)/sqrt(jac_x)/2;
        block_x{3,i} = block_x{3,i} + avg*jump';
        
        %Next edge integrals for functions who only share edge
        jump = -leg_edge_vals(:,1)/sqrt(jac_x);
        avg  = leg_edge_vals(:,2)/sqrt(jac_x)/2;
        block_x{4,i} = block_x{4,i} + avg*jump';
    end
    
    if i < num_x %upper edge integral
        %First edge integrals for functions who share same volume
        jump = leg_edge_vals(:,2)/sqrt(jac_x);
        avg  = leg_edge_vals(:,2)/sqrt(jac_x)/2;
        block_x{3,i} = block_x{3,i} + avg*jump';
        
        %Next edge integrals for functions who only share edge
        jump = leg_edge_vals(:,2)/sqrt(jac_x);
        avg  = leg_edge_vals(:,1)/sqrt(jac_x)/2;
        block_x{5,i} = block_x{5,i} + avg*jump';
    end
    
end

block_v = cell(10,num_v);
for i=1:num_v    
    %Create quadrature points on x cell
    quad_v = quad_ref*(v(i+1)-v(i))/2 + (v(i+1)+v(i))/2;
    
    %%%Create appropriate volume integrals in the v direction.    
    %(phi_i,phi_j)_{T_v}
    block_v{1,i} = (leg_vals/sqrt(jac_v))*(w_ref.*leg_vals'/sqrt(jac_v))*jac_v;
    %(\grad phi_i,phi_j)_{T_v}
    block_v{2,i} = (leg_vals/sqrt(jac_v))*(w_ref.*leg_der_vals'/(sqrt(jac_v)*jac_v))*jac_v;
    
    %Create edge integrals 
    block_v{3,i} = zeros(k+1);
    block_v{4,i} = zeros(k+1);
    block_v{5,i} = zeros(k+1);
    if i > 1 %lower edge integral
        %First edge integrals for functions who share same volume
        jump = -leg_edge_vals(:,1)/sqrt(jac_v);
        avg =   leg_edge_vals(:,1)/sqrt(jac_v)/2;
        block_v{3,i} = block_v{3,i} + avg*jump';
        
        %Next edge integrals for functions who only share edge
        jump = -leg_edge_vals(:,1)/sqrt(jac_v);
        avg  = leg_edge_vals(:,2)/sqrt(jac_v)/2;
        block_v{4,i} = block_v{4,i} + avg*jump';
    end
    
    if i < num_v %upper edge integral
        %First edge integrals for functions who share same volume
        jump = leg_edge_vals(:,2)/sqrt(jac_v);
        avg  = leg_edge_vals(:,2)/sqrt(jac_v)/2;
        block_v{3,i} = block_v{3,i} + avg*jump';
        
        %Next edge integrals for functions who only share edge
        jump = leg_edge_vals(:,2)/sqrt(jac_v);
        avg  = leg_edge_vals(:,1)/sqrt(jac_v)/2;
        block_v{5,i} = block_v{5,i} + avg*jump';
    end

    
end

%%Assemble stiffness matrix
% L = zeros((k+1)^2*(numel(x)-1)*(numel(v)-1));
% 
% for i=1:num_x
%     for j=1:num_v
%         blockstart = (k+1)^2*((i-1)*num_v+(j-1));
%         indices = blockstart+1:blockstart+(k+1)^2;
%         L(indices,indices) = L(indices,indices) + kron(block_x{1,i},block_v{1,j});
%         L(indices,indices) = L(indices,indices) + kron(block_x{2,i},block_v{2,j});
%         L(indices,indices) = L(indices,indices) + kron(block_x{3,i},block_v{3,j});
%         if j>1
%             L(indices-(k+1)^2,indices) = L(indices-(k+1)^2,indices) + kron(block_x{3,i},block_v{4,j});
%         end
%         if j<num_v
%             L(indices+(k+1)^2,indices) = L(indices+(k+1)^2,indices) + kron(block_x{3,i},block_v{5,j});
%         end
%     end
% end

%sparsify
I = zeros((k+1)^2*10*max([num_x num_v]),1);
J = zeros(size(I));
S = zeros(size(I));

count = 1;
dof = (k+1)^2*num_x*num_v;
for i=1:num_x
    for j=1:num_v
        blockstart = (k+1)^2*((i-1)*num_v+(j-1));
        indices = blockstart+1:blockstart+(k+1)^2;
        
        %All integrals involving tau_1 first
        temp = zeros((k+1)^2);
        %(\diff_x g,tau_1)_T
        temp = temp + kron( block_x{2,i},block_v{1,j});
        %<[g],{\tau_1}>_e (sharing volume)
        temp = temp - kron( block_x{3,i},block_v{1,j});
        J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
        I(count:count+(k+1)^4-1) = reshape(repmat(indices,1,(k+1)^2),[],1);
        S(count:count+(k+1)^4-1) = temp(:);
        count = count + (k+1)^4;
        %<[g],{\tau_1}>_e (sharing edge)
        if i>1
            indices_dn = indices - num_v*(k+1)^2;
            temp = zeros((k+1)^2);
            temp = temp - kron(block_x{4,i},block_v{1,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_dn,1,(k+1)^2),[],1);
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end
        if i<num_x
            indices_up = indices + num_v*(k+1)^2;
            temp = zeros((k+1)^2);
            temp = temp - kron(block_x{5,i},block_v{1,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_up,1,(k+1)^2),[],1);
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end
        
        %Now all integrals involving tau_2
        temp = zeros((k+1)^2);
        %(\diff_y g,tau_2)_T
        temp = temp + kron( block_x{1,i},block_v{2,j});
        %<[g],{\tau_1}>_e (sharing volume)
        temp = temp - kron( block_x{1,i},block_v{3,j});
        J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
        I(count:count+(k+1)^4-1) = reshape(repmat(indices,1,(k+1)^2),[],1)+dof;
        S(count:count+(k+1)^4-1) = temp(:);
        count = count + (k+1)^4;
        %<[g],{\tau_2}>_e (sharing edge)
        if j>1
            indices_dn = indices - (k+1)^2;
            temp = zeros((k+1)^2);
            temp = temp - kron(block_x{1,i},block_v{4,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_dn,1,(k+1)^2),[],1)+dof;
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end
        if j<num_v
            indices_up = indices + (k+1)^2;
            temp = zeros((k+1)^2);
            temp = temp - kron(block_x{1,i},block_v{5,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_up,1,(k+1)^2),[],1)+dof;
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end
    end
end
I(count+1:end) = [];
J(count+1:end) = [];
S(count+1:end) = [];

L = sparse(I,J,S,2*dof,dof);
end