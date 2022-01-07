function [L,Lx,Lv] = buildAdvectionMatrixConst(x,v,k,a)
%Stiffness matrix for a\cdot\grad u where a=(a_1,a_2) with a_1,a_2>0.
%Using Upwinding
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
        
        %Integrals sharing only edge
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

%sparsify
I = zeros((k+1)^2*10*max([num_x num_v]),1);
J = zeros(size(I));
S = zeros(size(I));

count = 1;
for i=1:num_x
    for j=1:num_v
        blockstart = (k+1)^2*((i-1)*num_v+(j-1));
        indices = blockstart+1:blockstart+(k+1)^2;
        
        %Volume integrals
        temp = zeros((k+1)^2);
        %-(a_1*u,diff_x v)_T
        temp = temp + kron(-a(1)*block_x{2,i},block_v{1,j});
        %-(a_2*u,diff_y v)_T
        temp = temp + kron(-a(2)*block_x{1,i},block_v{2,j});
        %Edge integrals that share volume
        %Integrals in x
        %Lower edge
        if i==1
            %Add nothing
        else
            temp = temp + kron(a(1)*block_x{4,i},block_v{1,j});
            temp = temp + kron(a(1)*block_x{9,i},block_v{1,j});
        end
        %Upper edge
        if i==num_x
            temp = temp + kron(a(1)*block_x{6,i},block_v{1,j});
        else
            temp = temp + kron(a(1)*block_x{6,i},block_v{1,j});
            temp = temp + kron(a(1)*block_x{11,i},block_v{1,j});
        end
        %Integrals in v
        %Lower edge
        if j==1
            %Add nothing
        else
            temp = temp + kron(a(2)*block_x{1,i},block_v{4,j});
            temp = temp + kron(a(2)*block_x{1,i},block_v{9,j});
        end
        %Upper edge
        if j==num_v
            temp = temp + kron(a(2)*block_x{1,i},block_v{6,j});
        else
            temp = temp + kron(a(2)*block_x{1,i},block_v{6,j});
            temp = temp + kron(a(2)*block_x{1,i},block_v{11,j});
        end
        J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
        I(count:count+(k+1)^4-1) = reshape(repmat(indices,1,(k+1)^2),[],1);
        S(count:count+(k+1)^4-1) = temp(:);
        count = count + (k+1)^4;
        %Edge integrals that only share edge
        if i>1
            indices_dn = indices - num_v*(k+1)^2;
            temp = zeros((k+1)^2);
            temp = temp + kron(a(1)*block_x{5,i},block_v{1,j});
            temp = temp + kron(a(1)*block_x{10,i},block_v{1,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_dn,1,(k+1)^2),[],1);
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end
        if i<num_x
            indices_up = indices + num_v*(k+1)^2;
            temp = zeros((k+1)^2);
            temp = temp + kron(a(1)*block_x{7,i},block_v{1,j});
            temp = temp + kron(a(1)*block_x{12,i},block_v{1,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_up,1,(k+1)^2),[],1);
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end
        if j>1
            indices_dn = indices - (k+1)^2;
            temp = zeros((k+1)^2);
            temp = temp + kron(a(2)*block_x{1,i},block_v{5,j});
            temp = temp + kron(a(2)*block_x{1,i},block_v{10,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_dn,1,(k+1)^2),[],1);
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end
        if j<num_v
            indices_up = indices + (k+1)^2;
            temp = zeros((k+1)^2);
            temp = temp + kron(a(2)*block_x{1,i},block_v{7,j});
            temp = temp + kron(a(2)*block_x{1,i},block_v{12,j});
            J(count:count+(k+1)^4-1) = repmat(indices,(k+1)^2,1);
            I(count:count+(k+1)^4-1) = reshape(repmat(indices_up,1,(k+1)^2),[],1);
            S(count:count+(k+1)^4-1) = temp(:);
            count = count + (k+1)^4;
        end
    end
end
if I(end) == 0
I(count:end) = [];
J(count:end) = [];
S(count:end) = [];
end

L = sparse(I,J,S);

%Create L_x
I1 = zeros(5*num_x*(k+1)^2,1);
J1 = I;
S1 = I;

I2 = I1;
J2 = J1;
S2 = S1;

count = [1,1];
for i=1:num_x
    blockstart = (k+1)*(i-1);
    indices = blockstart+1:blockstart+(k+1);
    
    %Grad
    grad = -block_x{2,i};
    if i==1
        %Add nothing
    else
        grad = grad + block_x{4,i};
        grad = grad + block_x{9,i};
    end
    %Upper edge
    if i==num_x
        grad = grad + block_x{6,i};
    else
        grad = grad + block_x{6,i};
        grad = grad + block_x{11,i};
    end    
    J1(count(1):count(1)+(k+1)^2-1) = repmat(indices,(k+1),1);
    I1(count(1):count(1)+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S1(count(1):count(1)+(k+1)^2-1) = grad(:);
    count(1) = count(1) + (k+1)^2;
    
    %Mass
    mass = block_x{1,i};
    J2(count(2):count(2)+(k+1)^2-1) = repmat(indices,(k+1),1);
    I2(count(2):count(2)+(k+1)^2-1) = reshape(repmat(indices,1,(k+1)),[],1);
    S2(count(2):count(2)+(k+1)^2-1) = mass(:);
    count(2) = count(2) + (k+1)^2;
    
    %More grad
    %Edge integrals that only share edge
    if i>1
        indices_dn = indices - (k+1);
        grad = block_x{5,i} + block_x{10,i};
        J1(count(1):count(1)+(k+1)^2-1) = repmat(indices,(k+1),1);
        I1(count(1):count(1)+(k+1)^2-1) = reshape(repmat(indices_dn,1,(k+1)),[],1);
        S1(count(1):count(1)+(k+1)^2-1) = grad(:);
        count(1) = count(1) + (k+1)^2;
    end
    if i<num_x
        indices_up = indices + (k+1);
        grad = block_x{7,i} + block_x{12,i};
        J1(count(1):count(1)+(k+1)^2-1) = repmat(indices,(k+1),1);
        I1(count(1):count(1)+(k+1)^2-1) = reshape(repmat(indices_up,1,(k+1)),[],1);
        S1(count(1):count(1)+(k+1)^2-1) = grad(:);
        count(1) = count(1) + (k+1)^2;
    end
end

I1(count(1):end) = [];
J1(count(1):end) = [];
S1(count(1):end) = [];
I2(count(2):end) = [];
J2(count(2):end) = [];
S2(count(2):end) = [];

Lx = {full(sparse(I1,J1,S1)),full(sparse(I2,J2,S2))};

end
