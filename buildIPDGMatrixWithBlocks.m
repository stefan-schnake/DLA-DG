function [Acell] = buildIPDGMatrixWithBlocks(x,v,k,kappa)
%Stiffness matrix Laplacian with homogeneous Dirichlet
% (\grad g , \grad z)_\W - <{\grad g},[z]>_{E_h} - <{\grad z},[g]>_{E_h} 
% + 1/h<[g],[z]>_{E_h}
%Each element T\in Th is T=T_x \times T_v.


num_x = numel(x)-1;
num_v = numel(v)-1;

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

[leg_vals,leg_der_vals,leg_edge_vals,leg_edge_der_vals] = buildLegendre(10,k);

block_x = cell(10,num_x);
for i=1:num_x   
    %Create quadrature points on x cell
    quad_x = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
    
    %%%Create appropriate integrals in the x direction.  
    %Each term in this decomposition is seperable so we can compute the
    %integrals in the x-direction and apply them to each v cell
    
    %(phi_i,phi_j)_{T_x}
    block_x{1,i} = (leg_vals/sqrt(jac_x))*(w_ref.*leg_vals'/sqrt(jac_x))*jac_x;
    %(\grad phi_i,\grad phi_j)_{T_x}
    block_x{2,i} = (leg_der_vals'/(jac_x*sqrt(jac_x)))*(w_ref.*leg_der_vals/(jac_x*sqrt(jac_x)))*jac_x;
    
    %Create edge integrals 
    block_x{3,i} = zeros(k+1);
    block_x{4,i} = zeros(k+1);
    block_x{5,i} = zeros(k+1);
    block_x{6,i} = zeros(k+1);
    block_x{7,i} = zeros(k+1);
    block_x{9,i} = zeros(k+1);
    block_x{10,i} = zeros(k+1);
    
    %%Lower edge integrals
    if i > 1
        %Integrals sharing volume
        %<{\grad_x phi_i},[phi_j]>
        I =  leg_edge_der_vals(:,1)/sqrt(jac_x^3)/2;
        J = -leg_edge_vals(:,1)/sqrt(jac_x);
        block_x{3,i} = block_x{3,i} + J*I';
        
        %<[phi_i],{\grad_x phi_j}>
        I = -leg_edge_vals(:,1)/sqrt(jac_x);
        J =  leg_edge_der_vals(:,1)/sqrt(jac_x^3)/2;
        block_x{3,i} = block_x{3,i} + J*I';
        
        %K/h<[phi_i],[phi_j]>
        I = -leg_edge_vals(:,1)/sqrt(jac_x);
        J = -leg_edge_vals(:,1)/sqrt(jac_x);
        block_x{4,i} = block_x{4,i} + kappa/(2*jac_x)*J*I';
        
        %Itegrals sharing only edge
        %<{\grad_x phi_i},[phi_j]>
        I = leg_edge_der_vals(:,1)/sqrt(jac_x^3)/2;
        J = leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{5,i} = block_x{5,i} + J*I';
        
        %<[phi_i],{\grad_x phi_j},>
        I = -leg_edge_vals(:,1)/sqrt(jac_x);
        J =  leg_edge_der_vals(:,2)/sqrt(jac_x^3)/2;
        block_x{5,i} = block_x{5,i} + J*I';
        
        %K/h<[phi_i],[phi_j]>
        I = -leg_edge_vals(:,1)/sqrt(jac_x);
        J =  leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{6,i} = block_x{6,i} + kappa/(2*jac_x)*J*I';
        
    else
        %Boundary integrals
        %Integrals sharing volume (negative is due to normal)
        %<\grad_x phi_i \cdot\nu,phi_j>
        I = -leg_edge_der_vals(:,1)/sqrt(jac_x^3);
        J =  leg_edge_vals(:,1)/sqrt(jac_x);
        block_x{3,i} = block_x{3,i} + J*I';
        
        %<phi_i,\grad_x phi_j \cdot\nu>
        I =  leg_edge_vals(:,1)/sqrt(jac_x);
        J = -leg_edge_der_vals(:,1)/sqrt(jac_x^3);
        block_x{3,i} = block_x{3,i} + J*I';
        
        %K/h<phi_i,phi_j>
        I = leg_edge_vals(:,1)/sqrt(jac_x);
        J = leg_edge_vals(:,1)/sqrt(jac_x);
        block_x{4,i} = block_x{4,i} + kappa/(2*jac_x)*J*I';
    end    
      
    %%Upper edge integrals  
    if i < num_x
        %Integrals sharing volume
        %<{\grad_x phi_i},[phi_j]>
        I = leg_edge_der_vals(:,2)/sqrt(jac_x^3)/2;
        J = leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{7,i} = block_x{7,i} + J*I';
        
        %<{\grad_x phi_i},[phi_j]>
        I = leg_edge_vals(:,2)/sqrt(jac_x);
        J = leg_edge_der_vals(:,2)/sqrt(jac_x^3)/2;
        block_x{7,i} = block_x{7,i} + J*I';
        
        %<K/h[phi_i],[phi_j]>
        I = leg_edge_vals(:,2)/sqrt(jac_x);
        J = leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{8,i} = block_x{8,i} + kappa/(2*jac_x)*J*I';
        
        %Integrals sharing only edge
        %<{\grad_x phi_i},[phi_j]>
        I = leg_edge_der_vals(:,2)/sqrt(jac_x^3)/2;
        J = -leg_edge_vals(:,1)/sqrt(jac_x);        
        block_x{9,i} = block_x{9,i} + J*I';
        
        %<[phi_i],{\grad_x phi_j}>
        I = leg_edge_vals(:,2)/sqrt(jac_x);        
        J = leg_edge_der_vals(:,1)/sqrt(jac_x^3)/2;
        block_x{9,i} = block_x{9,i} + J*I';
        
        %<K/h[phi_i],[phi_j]>
        I =  leg_edge_vals(:,2)/sqrt(jac_x);
        J = -leg_edge_vals(:,1)/sqrt(jac_x);        
        block_x{10,i} = block_x{12,i} + kappa/(2*jac_x)*J*I';
    else
        %Boundary integrals
        %Integrals sharing volume
        %<\grad_x phi_i\cdot\nu,phi_j>
        I = leg_edge_der_vals(:,2)/sqrt(jac_x^3);
        J = leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{7,i} = block_x{7,i} + J*I';
        
        %<phi_i,\grad_x phi_j\cdot\nu>
        I = leg_edge_vals(:,2)/sqrt(jac_x);
        J = leg_edge_der_vals(:,2)/sqrt(jac_x^3);
        block_x{7,i} = block_x{7,i} + J*I';
        
        %K/h<phi_i,phi_j>
        I = leg_edge_vals(:,2)/sqrt(jac_x);
        J = leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{8,i} = block_x{8,i} + kappa/(2*jac_x)*J*I';
    end

end

block_v = cell(10,num_v);
for i=1:num_v  
    %Create quadrature points on x cell
    quad_v = quad_ref*(v(i+1)-v(i))/2 + (v(i+1)+v(i))/2;
    
    %%%Create appropriate integrals in the v direction.  
    %Each term in this decomposition is seperable so we can compute the
    %integrals in the v-direction and apply them to each x cell
    
    %(phi_i,phi_j)_{T_v}
    block_v{1,i} = (leg_vals/sqrt(jac_v))*(w_ref.*leg_vals'/sqrt(jac_v))*jac_v;
    %(\grad phi_i,\grad phi_j)_{T_v}
    block_v{2,i} = (leg_der_vals'/(jac_v*sqrt(jac_v)))*(w_ref.*leg_der_vals/(jac_v*sqrt(jac_v)))*jac_v;
    
    %Create edge integrals 
    block_v{3,i} = zeros(k+1);
    block_v{4,i} = zeros(k+1);
    block_v{5,i} = zeros(k+1);
    block_v{6,i} = zeros(k+1);
    block_v{7,i} = zeros(k+1);
    block_v{9,i} = zeros(k+1);
    block_v{10,i} = zeros(k+1);
    
    %%Lower edge integrals
    if i > 1
        %Integrals sharing volume
        %<{\grad_v phi_i},[phi_j]>
        I =  leg_edge_der_vals(:,1)/sqrt(jac_v^3)/2;
        J = -leg_edge_vals(:,1)/sqrt(jac_v);
        block_v{3,i} = block_v{3,i} + J*I';
        
        %<[phi_i],{\grad_v phi_j}>
        I = -leg_edge_vals(:,1)/sqrt(jac_v);
        J =  leg_edge_der_vals(:,1)/sqrt(jac_v^3)/2;
        block_v{3,i} = block_v{3,i} + J*I';
        
        %K/h<[phi_i],[phi_j]>
        I = -leg_edge_vals(:,1)/sqrt(jac_v);
        J = -leg_edge_vals(:,1)/sqrt(jac_v);
        block_v{4,i} = block_v{4,i} + kappa/(2*jac_v)*J*I';
        
        %Itegrals sharing only edge
        %<{\grad_v phi_i},[phi_j]>
        I = leg_edge_der_vals(:,1)/sqrt(jac_v^3)/2;
        J = leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{5,i} = block_v{5,i} + J*I';
        
        %<[phi_i],{\grad_v phi_j},>
        I = -leg_edge_vals(:,1)/sqrt(jac_v);
        J =  leg_edge_der_vals(:,2)/sqrt(jac_v^3)/2;
        block_v{5,i} = block_v{5,i} + J*I';
        
        %K/h<[phi_i],[phi_j]>
        I = -leg_edge_vals(:,1)/sqrt(jac_v);
        J =  leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{6,i} = block_v{6,i} + kappa/(2*jac_v)*J*I';
        
    else
        %Boundary integrals
        %Integrals sharing volume (negative is due to normal)
        %<\grad_v phi_i \cdot\nu,phi_j>
        I = -leg_edge_der_vals(:,1)/sqrt(jac_v^3);
        J =  leg_edge_vals(:,1)/sqrt(jac_v);
        block_v{3,i} = block_v{3,i} + J*I';
        
        %<phi_i,\grad_v phi_j \cdot\nu>
        I =  leg_edge_vals(:,1)/sqrt(jac_v);
        J = -leg_edge_der_vals(:,1)/sqrt(jac_v^3);
        block_v{3,i} = block_v{3,i} + J*I';
        
        %K/h<phi_i,phi_j>
        I = leg_edge_vals(:,1)/sqrt(jac_v);
        J = leg_edge_vals(:,1)/sqrt(jac_v);
        block_v{4,i} = block_v{4,i} + kappa/(2*jac_v)*J*I';
    end    
      
    %%Upper edge integrals  
    if i < num_v
        %Integrals sharing volume
        %<{\grad_v phi_i},[phi_j]>
        I = leg_edge_der_vals(:,2)/sqrt(jac_v^3)/2;
        J = leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{7,i} = block_v{7,i} + J*I';
        
        %<{\grad_v phi_i},[phi_j]>
        I = leg_edge_vals(:,2)/sqrt(jac_v);
        J = leg_edge_der_vals(:,2)/sqrt(jac_v^3)/2;
        block_v{7,i} = block_v{7,i} + J*I';
        
        %<K/h[phi_i],[phi_j]>
        I = leg_edge_vals(:,2)/sqrt(jac_v);
        J = leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{8,i} = block_v{8,i} + kappa/(2*jac_v)*J*I';
        
        %Integrals sharing only edge
        %<{\grad_v phi_i},[phi_j]>
        I = leg_edge_der_vals(:,2)/sqrt(jac_v^3)/2;
        J = -leg_edge_vals(:,1)/sqrt(jac_v);        
        block_v{9,i} = block_v{9,i} + J*I';
        
        %<[phi_i],{\grad_v phi_j}>
        I = leg_edge_vals(:,2)/sqrt(jac_v);        
        J = leg_edge_der_vals(:,1)/sqrt(jac_v^3)/2;
        block_v{9,i} = block_v{9,i} + J*I';
        
        %<K/h[phi_i],[phi_j]>
        I =  leg_edge_vals(:,2)/sqrt(jac_v);
        J = -leg_edge_vals(:,1)/sqrt(jac_v);        
        block_v{10,i} = block_v{12,i} + kappa/(2*jac_v)*J*I';
    else
        %Boundary integrals
        %Integrals sharing volume
        %<\grad_v phi_i\cdot\nu,phi_j>
        I = leg_edge_der_vals(:,2)/sqrt(jac_v^3);
        J = leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{7,i} = block_v{7,i} + J*I';
        
        %<phi_i,\grad_v phi_j\cdot\nu>
        I = leg_edge_vals(:,2)/sqrt(jac_v);
        J = leg_edge_der_vals(:,2)/sqrt(jac_v^3);
        block_v{7,i} = block_v{7,i} + J*I';
        
        %K/h<phi_i,phi_j>
        I = leg_edge_vals(:,2)/sqrt(jac_v);
        J = leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{8,i} = block_v{8,i} + kappa/(2*jac_v)*J*I';
    end

end

clear I J

%%Term1 (pde in x times mass in v)
%Grad_X term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = block_x{2,i};
    grad = grad - block_x{3,i};
    grad = grad + block_x{4,i};
    grad = grad - block_x{7,i};
    grad = grad + block_x{8,i};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
    
    %More grad
    %Edge integrals that only share edge
    if i>1
        indices_dn = indices - (k+1);
        grad = -block_x{5,i} + block_x{6,i};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_dn,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
    if i<num_x
        indices_up = indices + (k+1);
        grad = -block_x{9,i} + block_x{10,i};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_up,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{1,1} = sparse(I,J,S,num_x*(k+1),num_x*(k+1));

%Mass_V term
count = 1;
I = zeros(5*num_v*(k+1)^2,1);
J = I;
S = I;
for j=1:num_v
    blockstart = (k+1)*(j-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    mass = block_v{1,j};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = mass(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{1,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));

%%Term5 (mass in x and pde in v)
%Mass_X term
count = 1;
I = zeros(5*num_v*(k+1)^2,1);
J = I;
S = I;
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    mass = block_x{1,i};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = mass(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{2,1} = sparse(I,J,S,num_x*(k+1),num_x*(k+1));

%Grad_V term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for j=1:num_v
    blockstart = (k+1)*(j-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = block_v{2,j};
    grad = grad - block_v{3,j};
    grad = grad + block_v{4,j};
    grad = grad - block_v{7,j};
    grad = grad + block_v{8,j};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
    
    
    %More grad
    %Edge integrals that only share edge
    if j>1
        indices_dn = indices - (k+1);
        grad = -block_v{5,j} + block_v{6,j};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_dn,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
    if j<num_x
        indices_up = indices + (k+1);
        grad = -block_v{9,j} + block_v{10,j};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_up,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{2,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));

end

