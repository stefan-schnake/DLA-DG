function [F] = evalBurgers2D(x,v,k,u)
%Evaluates weak form
%A(u;z) = -0.5*( (u,u)u,grad z )_\W + <\hat F(u),[z]>_E
%
%where \hat F(u) = {F(u)} + (max |u|)/2 [u].

num_x = numel(x)-1;
num_v = numel(v)-1;

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

[~, w_ref]  = lgwt(10,-1,1);

W = w_ref.*w_ref';

[leg_vals,leg_der_vals,leg_edge_vals,~] = buildLegendre(10,k);

F = zeros(num_x*num_v*(k+1)^2,1);

for i=1:num_x
    for j=1:num_v
        blockstart = (k+1)^2*((i-1)*num_v+(j-1));
        %%Volume integral       
        %Need u^2
        u_val = zeros(10);
        for ii=1:k+1
            for jj=1:k+1
                u_val = u_val + u(blockstart+(ii-1)*(k+1)+jj)*(leg_vals(ii,:)'*leg_vals(jj,:))/sqrt(jac_x*jac_v);
            end
        end
        for ii=1:k+1
            for jj=1:k+1
                %-(u^2/2,\grad_x z)_\W
                F(blockstart+(ii-1)*(k+1)+jj) = F(blockstart+(ii-1)*(k+1)+jj) - ...
                    .5*sum( u_val.^2.*W.*(leg_der_vals(ii,:)'*leg_vals(jj,:))/(sqrt(jac_x*jac_v)*jac_x) ,'all')*jac_x*jac_v;
                %-(u^2/2,\grad_v z)_\W
                F(blockstart+(ii-1)*(k+1)+jj) = F(blockstart+(ii-1)*(k+1)+jj) - ...
                    .5*sum( u_val.^2.*W.*(leg_vals(ii,:)'*leg_der_vals(jj,:))/(sqrt(jac_x*jac_v)*jac_v) ,'all')*jac_x*jac_v;
            end
        end
        
        u_max = norm(u_val(:),'inf');
        
        %%Edge integral
        %Bottom edge
        if j > 1
            u_edge = zeros(1,10);
            u_eg_mo = zeros(1,10); %Other edge
            u_vol_mo = zeros(10,10);
            mod = j-1;
            mod_bk = (k+1)^2*((i-1)*num_v+(mod-1));
            for ii=1:k+1
                for jj=1:k+1
                    u_edge  = u_edge  + u(blockstart+(ii-1)*(k+1)+jj)*leg_edge_vals(jj,1)*leg_vals(ii,:)/sqrt(jac_x*jac_v);
                    u_eg_mo = u_eg_mo + u(mod_bk    +(ii-1)*(k+1)+jj)*leg_edge_vals(jj,2)*leg_vals(ii,:)/sqrt(jac_x*jac_v);
                    u_vol_mo = u_vol_mo + u(mod_bk+(ii-1)*(k+1)+jj)*(leg_vals(ii,:)'*leg_vals(jj,:))/sqrt(jac_x*jac_v);
                end
            end        
            LF = max([norm(u_vol_mo(:),'inf'),u_max]);
            %{F(u)} + LF/2 [u]
            flux = 1/4*(u_edge.^2+u_eg_mo.^2) + LF/2*(-u_edge+u_eg_mo);   
            for ii=1:k+1
                for jj=1:k+1
                    jump = -leg_edge_vals(jj,1)*leg_vals(ii,:)/sqrt(jac_x*jac_v);
                    %<\hat F(u),[v]>
                    F(blockstart+(ii-1)*(k+1)+jj) = F(blockstart+(ii-1)*(k+1)+jj) + ...
                        flux*(jump'.*w_ref)*jac_x;
                end
            end
        end
        
        %Top edge
        if j < num_v
            u_edge = zeros(1,10);
            u_eg_mo = zeros(1,10); %Other edge
            u_vol_mo = zeros(10,10);
            mod = j+1;
            mod_bk = (k+1)^2*((i-1)*num_v+(mod-1));
            for ii=1:k+1
                for jj=1:k+1
                    u_edge  = u_edge  + u(blockstart+(ii-1)*(k+1)+jj)*leg_edge_vals(jj,2)*leg_vals(ii,:)/sqrt(jac_x*jac_v);
                    u_eg_mo = u_eg_mo + u(mod_bk    +(ii-1)*(k+1)+jj)*leg_edge_vals(jj,1)*leg_vals(ii,:)/sqrt(jac_x*jac_v);
                    u_vol_mo = u_vol_mo + u(mod_bk+(ii-1)*(k+1)+jj)*(leg_vals(ii,:)'*leg_vals(jj,:))/sqrt(jac_x*jac_v);
                end
            end
            LF = max([norm(u_vol_mo(:),'inf'),u_max]);
            %{F(u)} + LF/2 [u]
            flux = 1/4*(u_edge.^2+u_eg_mo.^2) + LF/2*(u_edge-u_eg_mo);
            for ii=1:k+1
                for jj=1:k+1
                    jump = leg_edge_vals(jj,2)*leg_vals(ii,:)/sqrt(jac_x*jac_v);
                    %<\hat F(u),[v]>
                    F(blockstart+(ii-1)*(k+1)+jj) = F(blockstart+(ii-1)*(k+1)+jj) + ...
                        flux*(jump'.*w_ref)*jac_x;
                end
            end
        else
            u_edge = zeros(1,10);
            for ii=1:k+1
                for jj=1:k+1
                    u_edge  = u_edge  + u(blockstart+(ii-1)*(k+1)+jj)*leg_edge_vals(jj,2)*leg_vals(ii,:)/sqrt(jac_x*jac_v);
                end
            end
            %{F(u)}
            flux = 1/2*(u_edge.^2);
            for ii=1:k+1
                for jj=1:k+1
                    jump = leg_edge_vals(jj,2)*leg_vals(ii,:)/sqrt(jac_x*jac_v);
                    %<\hat F(u),[v]>
                    F(blockstart+(ii-1)*(k+1)+jj) = F(blockstart+(ii-1)*(k+1)+jj) + ...
                        flux*(jump'.*w_ref)*jac_x;
                end
            end
        end
        
        %Left edge
        if i > 1
            u_edge = zeros(1,10);
            u_eg_mo = zeros(1,10); %Other edge
            u_vol_mo = zeros(10,10);
            mod = i-1;
            mod_bk = (k+1)^2*((mod-1)*num_v+(j-1));
            for ii=1:k+1
                for jj=1:k+1
                    u_edge  = u_edge  + u(blockstart+(ii-1)*(k+1)+jj)*leg_edge_vals(ii,1)*leg_vals(jj,:)/sqrt(jac_x*jac_v);
                    u_eg_mo = u_eg_mo + u(mod_bk    +(ii-1)*(k+1)+jj)*leg_edge_vals(ii,2)*leg_vals(jj,:)/sqrt(jac_x*jac_v);
                    u_vol_mo = u_vol_mo + u(mod_bk+(ii-1)*(k+1)+jj)*(leg_vals(ii,:)'*leg_vals(jj,:))/sqrt(jac_x*jac_v);
                end
            end        
            LF = max([norm(u_vol_mo(:),'inf'),u_max]);
            %{F(u)} + LF/2 [u]
            flux = 1/4*(u_edge.^2+u_eg_mo.^2) + LF/2*(-u_edge+u_eg_mo);     
            for ii=1:k+1
                for jj=1:k+1
                    jump = -leg_edge_vals(ii,1)*leg_vals(jj,:)/sqrt(jac_x*jac_v);
                    %<\hat F(u),[v]>
                    F(blockstart+(ii-1)*(k+1)+jj) = F(blockstart+(ii-1)*(k+1)+jj) + ...
                        flux*(jump'.*w_ref)*jac_v;
                end
            end
        end
        
        %Right edge
        if i < num_x
            u_edge = zeros(1,10);
            u_eg_mo = zeros(1,10); %Other edge
            u_vol_mo = zeros(10,10);
            mod = i+1; if mod > num_x; mod = 1; end
            mod_bk = (k+1)^2*((mod-1)*num_v+(j-1));
            for ii=1:k+1
                for jj=1:k+1
                    u_edge  = u_edge  + u(blockstart+(ii-1)*(k+1)+jj)*leg_edge_vals(ii,2)*leg_vals(jj,:)/sqrt(jac_x*jac_v);
                    u_eg_mo = u_eg_mo + u(mod_bk    +(ii-1)*(k+1)+jj)*leg_edge_vals(ii,1)*leg_vals(jj,:)/sqrt(jac_x*jac_v);
                    u_vol_mo = u_vol_mo + u(mod_bk+(ii-1)*(k+1)+jj)*(leg_vals(ii,:)'*leg_vals(jj,:))/sqrt(jac_x*jac_v);
                end
            end        
            LF = max([norm(u_vol_mo(:),'inf'),u_max]);
            %{F(u)} + LF/2 [u]
            flux = 1/4*(u_edge.^2+u_eg_mo.^2) + LF/2*(u_edge-u_eg_mo);     
            for ii=1:k+1
                for jj=1:k+1
                    jump = leg_edge_vals(ii,2)*leg_vals(jj,:)/sqrt(jac_x*jac_v);
                    %<\hat F(u),[v]>
                    F(blockstart+(ii-1)*(k+1)+jj) = F(blockstart+(ii-1)*(k+1)+jj) + ...
                        flux*(jump'.*w_ref)*jac_v;
                end
            end
        else
            u_edge = zeros(1,10);
            for ii=1:k+1
                for jj=1:k+1
                    u_edge  = u_edge  + u(blockstart+(ii-1)*(k+1)+jj)*leg_edge_vals(ii,2)*leg_vals(jj,:)/sqrt(jac_x*jac_v);
                end
            end        
            %{F(u)}
            flux = 1/2*(u_edge.^2);   
            for ii=1:k+1
                for jj=1:k+1
                    jump = leg_edge_vals(ii,2)*leg_vals(jj,:)/sqrt(jac_x*jac_v);
                    %<\hat F(u),[v]>
                    F(blockstart+(ii-1)*(k+1)+jj) = F(blockstart+(ii-1)*(k+1)+jj) + ...
                        flux*(jump'.*w_ref)*jac_v;
                end
            end
        end
        
        
        
    end
end

end
