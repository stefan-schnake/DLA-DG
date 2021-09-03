function [leg_vals,leg_der_vals,leg_edge_vals,leg_der_edge_vals] = buildLegendre(N,k)
%Create legendre polynomial quadrature points 
[quad_ref, ~]  = lgwt(N,-1,1); 
quad_ref = [-1;quad_ref;1]';

%Get legendre polynomial values at quad_ref by recursion
leg_val = zeros(k+1,numel(quad_ref));
leg_val(1,:) = 0.*quad_ref+1;
leg_val(2,:) = quad_ref;
for i=2:k
    n = i-1;
    leg_val(i+1,:) = ( (2*n+1)*quad_ref.*leg_val(i,:)-n*leg_val(i-1,:) )/(n+1);
end

leg_der_val = zeros(k+1,numel(quad_ref));
leg_der_val(1,:) = 0*quad_ref;
leg_der_val(2,:) = 0*quad_ref+1;
for i=2:k
    n = i-1;
    leg_der_val(i+1,:) = (n+1)*leg_val(i,:) + quad_ref.*leg_der_val(i,:);
end


%Normalize in L^2(-1,1) and split
leg_vals = zeros(k+1,numel(quad_ref)-2);
leg_der_vals = zeros(k+1,numel(quad_ref)-2);
leg_edge_vals = zeros(k+1,2);
leg_der_edge_vals = zeros(k+1,2);
for i=1:k+1
    leg_val(i,:) = leg_val(i,:)*sqrt((2*(i-1)+1)/2);
    leg_der_val(i,:) = leg_der_val(i,:)*sqrt((2*(i-1)+1)/2);
    leg_vals(i,:) = leg_val(i,2:end-1);
    leg_der_vals(i,:) = leg_der_val(i,2:end-1);
    leg_edge_vals(i,:) = leg_val(i,[1,numel(quad_ref)]);
    leg_der_edge_vals(i,:) = leg_der_val(i,[1,numel(quad_ref)]);
end


end

