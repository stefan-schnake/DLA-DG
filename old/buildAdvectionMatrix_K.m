function [F] = buildAdvectionMatrix_K(x,v,k,K,D)
%DLA update on K for a\cdot\grad u where a=(-y,x)

num_x = numel(x)-1;
num_v = numel(v)-1;

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

[leg_vals,leg_der_vals,leg_edge_vals,~] = buildLegendre(10,k);

F = K*0;
r = size(K,2);

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



for i=1:num_x
    indx = (i-1)*(k+1)+1:i*(k+1);
    for j=1:num_v
        indv = (j-1)*(k+1)+1:j*(k+1);
        %Volume integrals
        for p=1:r
            for q=1:r
                F(indx,q) = F(indx,q) - block_x{2,i}*K(indx,p)*(D(indv,q)'*block_v{3,j}*D(indv,p));
                F(indx,q) = F(indx,q) - block_x{3,i}*K(indx,p)*(D(indv,q)'*block_v{2,j}*D(indv,p));
            end
        end
        %Edge integrals that share volume
        %Lower edge
        if i==1
            if (v(j+1)+v(j)) > 0 %Outflow
                for p=1:r
                    for q=1:r
                        F(indx,q) = F(indx,q) + block_x{4,i}*K(indx,p)*(D(indv,q)'*block_v{3,j}*D(indv,p));
                    end
                end 
            end
        else
            for p=1:r
                for q=1:r
                    F(indx,q) = F(indx,q) + block_x{4,i}*K(indx,p)*(D(indv,q)'*block_v{3,j}*D(indv,p));
                    F(indx,q) = F(indx,q) + block_x{9,i}*K(indx,p)*(D(indv,q)'*block_v{8,j}*D(indv,p));
                end
            end 
        end
        %Upper edge
        if i==num_x
            if (v(j+1)+v(j)) < 0 %Outflow
                for p=1:r
                    for q=1:r
                        F(indx,q) = F(indx,q) + block_x{6,i}*K(indx,p)*(D(indv,q)'*block_v{3,j}*D(indv,p));
                    end
                end 
            end
        else
            for p=1:r
                for q=1:r
                    F(indx,q) = F(indx,q) + block_x{6,i}*K(indx,p)*(D(indv,q)'*block_v{3,j}*D(indv,p));
                    F(indx,q) = F(indx,q) + block_x{11,i}*K(indx,p)*(D(indv,q)'*block_v{8,j}*D(indv,p));
                end
            end
        end
        %Integrals in v
        %Lower edge
        if j==1
            if (x(i+1)+x(i)) < 0 %Outflow
                for p=1:r
                    for q=1:r
                        F(indx,q) = F(indx,q) + block_x{3,i}*K(indx,p)*(D(indv,q)'*block_v{4,j}*D(indv,p));
                    end
                end 
            end
        else
            for p=1:r
                for q=1:r
                    F(indx,q) = F(indx,q) + block_x{3,i}*K(indx,p)*(D(indv,q)'*block_v{4,j}*D(indv,p));
                    F(indx,q) = F(indx,q) + block_x{8,i}*K(indx,p)*(D(indv,q)'*block_v{9,j}*D(indv,p));
                end
            end
        end
        %Upper edge
        if j==num_v
            if (x(i+1)+x(i)) > 0 %Outflow
                for p=1:r
                    for q=1:r
                        F(indx,q) = F(indx,q) + block_x{3,i}*K(indx,p)*(D(indv,q)'*block_v{6,j}*D(indv,p));
                    end
                end 
            end
        else
            for p=1:r
                for q=1:r
                    F(indx,q) = F(indx,q) + block_x{3,i}*K(indx,p)*(D(indv,q)'*block_v{6,j}*D(indv,p));
                    F(indx,q) = F(indx,q) + block_x{8,i}*K(indx,p)*(D(indv,q)'*block_v{11,j}*D(indv,p));
                end
            end
        end
        
        %Edge integrals that only share edge
        %Left edge
        if i > 1
            indxmod = indx - (k+1);
            for p=1:r
                for q=1:r
                    F(indxmod,q) = F(indxmod,q) + block_x{5,i}*K(indx,p)*(D(indv,q)'*block_v{3,j}*D(indv,p));
                    F(indxmod,q) = F(indxmod,q) + block_x{10,i}*K(indx,p)*(D(indv,q)'*block_v{8,j}*D(indv,p));
                end
            end
        end
        %Right edge
        if i < num_x
            indxmod = indx + (k+1);
            for p=1:r
                for q=1:r
                    F(indxmod,q) = F(indxmod,q) + block_x{7,i}*K(indx,p)*(D(indv,q)'*block_v{3,j}*D(indv,p));
                    F(indxmod,q) = F(indxmod,q) + block_x{12,i}*K(indx,p)*(D(indv,q)'*block_v{8,j}*D(indv,p));
                end
            end
        end
        %Bottom edge
        if j > 1
            indvmod = indv - (k+1);
            for p=1:r
                for q=1:r
                    F(indx,q) = F(indx,q) + block_x{3,i}*K(indx,p)*(D(indvmod,q)'*block_v{5,j}*D(indv,p));
                    F(indx,q) = F(indx,q) + block_x{8,i}*K(indx,p)*(D(indvmod,q)'*block_v{10,j}*D(indv,p));
                end
            end
        end
        %Top edge
        if j < num_v
            indvmod = indv + (k+1);
            for p=1:r
                for q=1:r
                    F(indx,q) = F(indx,q) + block_x{3,i}*K(indx,p)*(D(indvmod,q)'*block_v{7,j}*D(indv,p));
                    F(indx,q) = F(indx,q) + block_x{8,i}*K(indx,p)*(D(indvmod,q)'*block_v{12,j}*D(indv,p));
                end
            end
        end
    end
end

end