function [Acell] = buildConstantAdvectionMatrixWithBlocks(x,v,k,a)
%Stiffness matrix for a\cdot\grad u where a is given and BCs are periodic
%Using average gradients for now
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
    %(phi_i,\grad phi_j)_{T_x}
    block_x{2,i} = (leg_der_vals/(jac_x*sqrt(jac_x)))*(w_ref.*leg_vals'/sqrt(jac_x))*jac_x;
    
    %Create edge integrals 
    block_x{4,i} = zeros(k+1);
    block_x{5,i} = zeros(k+1);
    block_x{6,i} = zeros(k+1);
    block_x{7,i} = zeros(k+1);
    block_x{8,i} = zeros(k+1);
    block_x{9,i} = zeros(k+1);
    
    %%Lower edge integrals
    
    %Integrals sharing volume
    %<{phi_i},[phi_j]>
    avg =   leg_edge_vals(:,1)/sqrt(jac_x)/2;
    jump = -leg_edge_vals(:,1)/sqrt(jac_x);
    block_x{4,i} = block_x{4,i} + jump*avg';

    %<[phi_i],[phi_j]>
    avg =  -leg_edge_vals(:,1)/sqrt(jac_x);
    jump = -leg_edge_vals(:,1)/sqrt(jac_x);
    block_x{5,i} = block_x{5,i} + jump*avg';

    %Itegrals sharing only edge
    %<{phi_i},[phi_j]>
    avg  = leg_edge_vals(:,1)/sqrt(jac_x)/2;
    jump = leg_edge_vals(:,2)/sqrt(jac_x);
    block_x{6,i} = block_x{6,i} + jump*avg';

    %<[phi_i],[phi_j]>
    avg =  -leg_edge_vals(:,1)/sqrt(jac_x);
    jump =  leg_edge_vals(:,2)/sqrt(jac_x);
    block_x{7,i} = block_x{7,i} + jump*avg';

      
    %%Upper edge integrals  

    %Integrals sharing volume
    %<{phi_i},[phi_j]>
    avg  = leg_edge_vals(:,2)/sqrt(jac_x)/2;
    jump = leg_edge_vals(:,2)/sqrt(jac_x);
    block_x{4,i} = block_x{4,i} + jump*avg';

    %<[phi_i],[phi_j]>
    avg  = leg_edge_vals(:,2)/sqrt(jac_x);
    jump = leg_edge_vals(:,2)/sqrt(jac_x);
    block_x{5,i} = block_x{5,i} + jump*avg';

    %Integrals sharing only edge
    %<{phi_i},[phi_j]>
    avg  = leg_edge_vals(:,2)/sqrt(jac_x)/2;
    jump = -leg_edge_vals(:,1)/sqrt(jac_x);        
    block_x{8,i} = block_x{8,i} + jump*avg';

    %<[phi_i],[phi_j]>
    avg  = leg_edge_vals(:,2)/sqrt(jac_x);
    jump = -leg_edge_vals(:,1)/sqrt(jac_x);        
    block_x{9,i} = block_x{9,i} + jump*avg';

end

block_v = cell(10,num_v);
for i=1:num_v    
    %Create quadrature points on x cell
    quad_v = quad_ref*(v(i+1)-v(i))/2 + (v(i+1)+v(i))/2;
    
    %%%Create appropriate volume integrals in the v direction.    
    %(phi_i,phi_j)_{T_v}
    block_v{1,i} = (leg_vals/sqrt(jac_v))*(w_ref.*leg_vals'/sqrt(jac_v))*jac_v;
    %(phi_i,\grad phi_j)_{T_v}
    block_v{2,i} = (leg_der_vals/(sqrt(jac_v)*jac_v))*(w_ref.*leg_vals'/sqrt(jac_v))*jac_v;
    
    %Create edge integrals 
    block_v{4,i} = zeros(k+1);
    block_v{5,i} = zeros(k+1);
    block_v{6,i} = zeros(k+1);
    block_v{7,i} = zeros(k+1);
    block_v{8,i} = zeros(k+1);
    block_v{9,i} = zeros(k+1);
    
    %%Lower edge integrals
    
    %Integrals sharing volume
    %<{phi_i},[phi_j]>
    avg =   leg_edge_vals(:,1)/sqrt(jac_v)/2;
    jump = -leg_edge_vals(:,1)/sqrt(jac_v);
    block_v{4,i} = block_v{4,i} + jump*avg';

    %<[phi_i],[phi_j]>
    avg =  -leg_edge_vals(:,1)/sqrt(jac_v);
    jump = -leg_edge_vals(:,1)/sqrt(jac_v);
    block_v{5,i} = block_v{5,i} + jump*avg';

    %Itegrals sharing only edge
    %<{phi_i},[phi_j]>
    avg  = leg_edge_vals(:,1)/sqrt(jac_v)/2;
    jump = leg_edge_vals(:,2)/sqrt(jac_v);
    block_v{6,i} = block_v{6,i} + jump*avg';

    %<[phi_i],[phi_j]>
    avg  = -leg_edge_vals(:,1)/sqrt(jac_v);
    jump = leg_edge_vals(:,2)/sqrt(jac_v);
    block_v{7,i} = block_v{7,i} + jump*avg';
      
    %%Upper edge integrals
    
    %Integrals sharing volume
    %<{phi_i},[phi_j]>
    avg  = leg_edge_vals(:,2)/sqrt(jac_v)/2;
    jump = leg_edge_vals(:,2)/sqrt(jac_v);
    block_v{4,i} = block_v{4,i} + jump*avg';

    %<[phi_i],[phi_j]>
    avg  = leg_edge_vals(:,2)/sqrt(jac_v);
    jump = leg_edge_vals(:,2)/sqrt(jac_v);
    block_v{5,i} = block_v{5,i} + jump*avg';

    %Integrals sharing only edge
    %<{phi_i},[phi_j]>
    avg  = leg_edge_vals(:,2)/sqrt(jac_v)/2;
    jump = -leg_edge_vals(:,1)/sqrt(jac_v);        
    block_v{8,i} = block_v{8,i} + jump*avg';

    %<[phi_i],[phi_j]>
    avg  = leg_edge_vals(:,2)/sqrt(jac_v);
    jump = -leg_edge_vals(:,1)/sqrt(jac_v);        
    block_v{9,i} = block_v{9,i} + jump*avg';
    
end


%%Term1 (PDE down in x * MASS in v)
%Grad_X term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = -a(1)*block_x{2,i} + a(1)*block_x{4,i} + abs(a(1))/2*block_x{5,i};
    
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
    
    
    %More grad
    %Edge integrals that only share edge
    
    %Lower edge
    if i>1
        indices_dn = indices - (k+1);
    else
        indices_dn = indices - (k+1) + num_x*(k+1);
    end
    grad = a(1)*block_x{6,i} + abs(a(1))/2*block_x{7,i};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices_dn,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
    if i<num_x
        indices_up = indices + (k+1);
    else
        indices_up = indices + (k+1) - num_x*(k+1);
    end        
    grad = a(1)*block_x{8,i} + abs(a(1))/2*block_x{9,i};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices_up,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;

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


%%Term2 (MASS in x and PDE upwind in v)
%Mass_X term
count = 1;
I = zeros(5*num_v*(k+1)^2,1);
J = I;
S = I;
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Mass
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
    grad = -a(2)*block_v{2,j} + a(2)*block_v{4,j} + abs(a(2))/2*block_v{5,j};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
    
    
    %More grad
    %Edge integrals that only share edge
    if j>1
        indices_dn = indices - (k+1);
    else
        indices_dn = indices - (k+1) + num_v*(k+1);
    end        
    grad = a(2)*block_v{6,j} + abs(a(2))/2*block_v{7,j};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices_dn,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;

    if j<num_x
        indices_up = indices + (k+1);
    else
        indices_up = indices + (k+1) - num_x*(k+1);
    end
    grad = a(2)*block_v{8,j} + abs(a(2))/2*block_v{9,j};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices_up,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{2,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));


end
