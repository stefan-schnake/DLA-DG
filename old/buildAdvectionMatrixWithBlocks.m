function [Acell] = buildAdvectionMatrixWithBlocks(x,v,k)
%Stiffness matrix for a\cdot\grad u where a=(-y,x) with 
%appropriate BCs.
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
    %(phi_i,x*phi_j)_{T_x}
    block_x{3,i} = (leg_vals/sqrt(jac_x).*quad_x)*(w_ref.*leg_vals'/sqrt(jac_x))*jac_x;
    %(1/2)*(phi_i,|x|*phi_j}_{T_x}
    block_x{8,i} = 1/2*(leg_vals/sqrt(jac_x).*abs(quad_x))*(w_ref.*leg_vals'/sqrt(jac_x))*jac_x;
    
    %Create edge integrals 
    block_x{4,i} = zeros(k+1);
    block_x{5,i} = zeros(k+1);
    block_x{6,i} = zeros(k+1);
    block_x{7,i} = zeros(k+1);
    block_x{9,i} = zeros(k+1);
    block_x{10,i} = zeros(k+1);
    block_x{11,i} = zeros(k+1);
    block_x{12,i} = zeros(k+1);
    
    %%Lower edge integrals
    if i > 1
        %Integrals sharing volume
        %<{phi_i},[phi_j]>
        avg =   leg_edge_vals(:,1)/sqrt(jac_x)/2;
        jump = -leg_edge_vals(:,1)/sqrt(jac_x);
        block_x{4,i} = block_x{4,i} + jump*avg';
        
        %<[phi_i],[phi_j]>
        avg =  -leg_edge_vals(:,1)/sqrt(jac_x);
        jump = -leg_edge_vals(:,1)/sqrt(jac_x);
        block_x{9,i} = block_x{9,i} + jump*avg';
        
        %Itegrals sharing only edge
        %<{phi_i},[phi_j]>
        avg  = leg_edge_vals(:,1)/sqrt(jac_x)/2;
        jump = leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{5,i} = block_x{5,i} + jump*avg';
        
        %<[phi_i],[phi_j]>
        avg =  -leg_edge_vals(:,1)/sqrt(jac_x);
        jump =  leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{10,i} = block_x{10,i} + jump*avg';
        
    else
        %Integrals sharing volume (negative is due to normal)
        %<phi_i,phi_j>_{E^-}
        jump = leg_edge_vals(:,1)/sqrt(jac_x);
        avg  = leg_edge_vals(:,1)/sqrt(jac_x);
        block_x{4,i} = block_x{4,i} - jump*avg';
    end    
      
    %%Upper edge integrals  
    if i < num_x
        %Integrals sharing volume
        %<{phi_i},[phi_j]>
        avg  = leg_edge_vals(:,2)/sqrt(jac_x)/2;
        jump = leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{6,i} = block_x{6,i} + jump*avg';
        
        %<[phi_i],[phi_j]>
        avg  = leg_edge_vals(:,2)/sqrt(jac_x);
        jump = leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{11,i} = block_x{11,i} + jump*avg';
        
        %Integrals sharing only edge
        %<{phi_i},[phi_j]>
        avg  = leg_edge_vals(:,2)/sqrt(jac_x)/2;
        jump = -leg_edge_vals(:,1)/sqrt(jac_x);        
        block_x{7,i} = block_x{7,i} + jump*avg';
        
        %<[phi_i],[phi_j]>
        avg  = leg_edge_vals(:,2)/sqrt(jac_x);
        jump = -leg_edge_vals(:,1)/sqrt(jac_x);        
        block_x{12,i} = block_x{12,i} + jump*avg';
    else
        %Integrals sharing volume
        %<phi_i,phi_j>_{E^+}
        jump = leg_edge_vals(:,2)/sqrt(jac_x);
        avg  = leg_edge_vals(:,2)/sqrt(jac_x);
        block_x{6,i} = block_x{6,i} + jump*avg';
    end

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
    %(phi_i*(-y),phi_j)_{T_v}
    block_v{3,i} = -(leg_vals/sqrt(jac_v).*quad_v)*(w_ref.*leg_vals'/sqrt(jac_v))*jac_v;
    %(1/2)*(phi_i*|y|,phi_j)_{T_v}
    block_v{8,i} = 1/2*(leg_vals/sqrt(jac_v).*abs(quad_v))*(w_ref.*leg_vals'/sqrt(jac_v))*jac_v;
    
    
    %Create edge integrals 
    block_v{4,i} = zeros(k+1);
    block_v{5,i} = zeros(k+1);
    block_v{6,i} = zeros(k+1);
    block_v{7,i} = zeros(k+1);
    block_v{9,i} = zeros(k+1);
    block_v{10,i} = zeros(k+1);
    block_v{11,i} = zeros(k+1);
    block_v{12,i} = zeros(k+1);
    
    %%Lower edge integrals
    if i > 1
        %Integrals sharing volume
        %<{phi_i},[phi_j]>
        avg =   leg_edge_vals(:,1)/sqrt(jac_v)/2;
        jump = -leg_edge_vals(:,1)/sqrt(jac_v);
        block_v{4,i} = block_v{4,i} + jump*avg';
        
        %<[phi_i],[phi_j]>
        avg =  -leg_edge_vals(:,1)/sqrt(jac_v);
        jump = -leg_edge_vals(:,1)/sqrt(jac_v);
        block_v{9,i} = block_v{9,i} + jump*avg';
        
        %Itegrals sharing only edge
        %<{phi_i},[phi_j]>
        avg  = leg_edge_vals(:,1)/sqrt(jac_v)/2;
        jump = leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{5,i} = block_v{5,i} + jump*avg';
        
        %<[phi_i],[phi_j]>
        avg  = -leg_edge_vals(:,1)/sqrt(jac_v);
        jump = leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{10,i} = block_v{10,i} + jump*avg';
    else
        %Integrals sharing volume (negative is due to normal)
        %<phi_i,phi_j>_{E^-}
        jump = leg_edge_vals(:,1)/sqrt(jac_v);
        avg  = leg_edge_vals(:,1)/sqrt(jac_v);
        block_v{4,i} = block_v{4,i} - jump*avg';
    end    
      
    %%Upper edge integrals
    
    if i < num_v
        %Integrals sharing volume
        %<{phi_i},[phi_j]>
        avg  = leg_edge_vals(:,2)/sqrt(jac_v)/2;
        jump = leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{6,i} = block_v{6,i} + jump*avg';
        
        %<[phi_i],[phi_j]>
        avg  = leg_edge_vals(:,2)/sqrt(jac_v);
        jump = leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{11,i} = block_v{11,i} + jump*avg';
        
        %Integrals sharing only edge
        %<{phi_i},[phi_j]>
        avg  = leg_edge_vals(:,2)/sqrt(jac_v)/2;
        jump = -leg_edge_vals(:,1)/sqrt(jac_v);        
        block_v{7,i} = block_v{7,i} + jump*avg';
        
        %<[phi_i],[phi_j]>
        avg  = leg_edge_vals(:,2)/sqrt(jac_v);
        jump = -leg_edge_vals(:,1)/sqrt(jac_v);        
        block_v{12,i} = block_v{12,i} + jump*avg';
    else
        %Integrals sharing volume
        %<phi_i,phi_j>_{E^+}
        jump = leg_edge_vals(:,2)/sqrt(jac_v);
        avg  = leg_edge_vals(:,2)/sqrt(jac_v);
        block_v{6,i} = block_v{6,i} + jump*avg';
    end
    
end


%%Term1 (PDE in x on volume and interior edges)
%Grad_X term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = -block_x{2,i};
    %Lower Edge
    if i==1
        %Add nothing
    else
        grad = grad + block_x{4,i};
    end
    %Upper edge
    if i==num_x
        %Add nothing
    else
        grad = grad + block_x{6,i};
    end    
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
    
    
    %More grad
    %Edge integrals that only share edge
    if i>1
        indices_dn = indices - (k+1);
        grad = block_x{5,i};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_dn,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
    if i<num_x
        indices_up = indices + (k+1);
        grad = block_x{7,i};
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
    mass = block_v{3,j};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = mass(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{1,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));

%%Term2 (penalty terms in x)
%Grad_X term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = zeros(k+1);
    %Lower Edge
    if i==1
        %Add nothing
    else
        grad = grad + block_x{9,i};
    end
    %Upper edge
    if i==num_x
        %Add nothing
    else
        grad = grad + block_x{11,i};
    end    
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
    
    
    %More grad
    %Edge integrals that only share edge
    if i>1
        indices_dn = indices - (k+1);
        grad = block_x{10,i};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_dn,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
    if i<num_x
        indices_up = indices + (k+1);
        grad = block_x{12,i};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_up,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{2,1} = sparse(I,J,S,num_x*(k+1),num_x*(k+1));

%Mass_V term
count = 1;
I = zeros(5*num_v*(k+1)^2,1);
J = I;
S = I;
for j=1:num_v
    blockstart = (k+1)*(j-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    mass = block_v{8,j};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = mass(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{2,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));

%%Term3 (pde outflow in x on bottom)
%Grad_X term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for i=1:1
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = block_x{4,i}; 
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{3,1} = sparse(I,J,S,num_x*(k+1),num_x*(k+1));

%Mass_V term
count = 1;
I = zeros(5*num_v*(k+1)^2,1);
J = I;
S = I;
for j=1:num_v
    blockstart = (k+1)*(j-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    if (v(j+1)+v(j)) > 0 %Outflow
        mass = block_v{3,j};
    else
        mass = zeros(k+1);
    end
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = mass(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{3,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));

%%Term4 (pde outflow in x on top)
%Grad_X term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for i=num_x:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = block_x{6,i}; 
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{4,1} = sparse(I,J,S,num_x*(k+1),num_x*(k+1));

%Mass_V term
count = 1;
I = zeros(5*num_v*(k+1)^2,1);
J = I;
S = I;
for j=1:num_v
    blockstart = (k+1)*(j-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %mass
    if (v(j+1)+v(j)) < 0 %Outflow
        mass = block_v{3,j};
    else
        mass = zeros(k+1);
    end
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = mass(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{4,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));

%%Term5 (PDE in v on volume and interior edges)
%Mass_X term
count = 1;
I = zeros(5*num_v*(k+1)^2,1);
J = I;
S = I;
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    mass = block_x{3,i};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = mass(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{5,1} = sparse(I,J,S,num_x*(k+1),num_x*(k+1));

%Grad_V term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for j=1:num_v
    blockstart = (k+1)*(j-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = -block_v{2,j};
    %Lower Edge
    if j==1
        %Add nothing
    else
        grad = grad + block_v{4,j};
    end
    %Upper edge
    if j==num_v
        %Add nothing
    else
        grad = grad + block_v{6,j};
    end    
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
    
    
    %More grad
    %Edge integrals that only share edge
    if j>1
        indices_dn = indices - (k+1);
        grad = block_v{5,j};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_dn,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
    if j<num_x
        indices_up = indices + (k+1);
        grad = block_v{7,j};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_up,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{5,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));

%%Term6 (jump terms in v on interior edges)
%Mass_X term
count = 1;
I = zeros(5*num_v*(k+1)^2,1);
J = I;
S = I;
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    mass = block_x{8,i};
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = mass(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{6,1} = sparse(I,J,S,num_x*(k+1),num_x*(k+1));

%Grad_V term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for j=1:num_v
    blockstart = (k+1)*(j-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = zeros(k+1);
    %Lower Edge
    if j==1
        %Add nothing
    else
        grad = grad + block_v{9,j};
    end
    %Upper edge
    if j==num_v
        %Add nothing
    else
        grad = grad + block_v{11,j};
    end    
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
    
    
    %More grad
    %Edge integrals that only share edge
    if j>1
        indices_dn = indices - (k+1);
        grad = block_v{10,j};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_dn,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
    if j<num_x
        indices_up = indices + (k+1);
        grad = block_v{12,j};
        J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
        I(count:count+(k+1)^2-1) = reshape(repmat(indices_up,1,(k+1)),[],1);
        S(count:count+(k+1)^2-1) = grad(:);
        count = count + (k+1)^2;
    end
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{6,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));

%%Term7 (outflow in x on left edge)
%Mass_X term
count = 1;
I = zeros(5*num_v*(k+1)^2,1);
J = I;
S = I;
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    if (x(i+1)+x(i)) < 0
        mass = block_x{3,i};
    else
        mass = zeros(k+1);
    end
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = mass(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{7,1} = sparse(I,J,S,num_x*(k+1),num_x*(k+1));

%Grad_V term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for j=1:1
    blockstart = (k+1)*(j-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = block_v{4,j};  
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{7,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));

%%Term8 (outflow in x on left edge)
%Mass_X term
count = 1;
I = zeros(5*num_v*(k+1)^2,1);
J = I;
S = I;
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    if (x(i+1)+x(i)) > 0
        mass = block_x{3,i};
    else
        mass = zeros(k+1);
    end
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = mass(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{8,1} = sparse(I,J,S,num_x*(k+1),num_x*(k+1));

%Grad_V term
count = 1;
I = zeros(5*num_x*(k+1)^2,1);
J = I;
S = I;
for j=num_v:num_v
    blockstart = (k+1)*(j-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = block_v{6,j};  
    J(count:count+(k+1)^2-1) = repmat(indices,(k+1),1);
    I(count:count+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S(count:count+(k+1)^2-1) = grad(:);
    count = count + (k+1)^2;
end
I(count:end) = []; J(count:end) = []; S(count:end) = [];
Acell{8,2} = sparse(I,J,S,num_v*(k+1),num_v*(k+1));


end
