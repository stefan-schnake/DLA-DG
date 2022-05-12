function [] = plotVecDualDiscrete(x,v,k,u,u2,myhist)
%Plots u


rel = false;


num_x = numel(x)-1;
num_v = numel(v)-1;
poi = 1;

label_x = 'x';
label_y = 'v'; 

[quad_ref, w_ref]  = lgwt(poi,-1,1);
quad_ref = quad_ref';
quad_ref = fliplr(quad_ref);

[leg_vals,~,~,~] = buildLegendre(poi,k);
leg_vals = fliplr(leg_vals);

jac_x = (x(2)-x(1))/2;
jac_v = (v(2)-v(1))/2;

X = zeros(poi*num_x,1);
V = zeros(poi*num_v,1);
for i=1:num_x
    X((i-1)*poi+1:i*poi) = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
end
for i=1:num_v
    V((i-1)*poi+1:i*poi) = quad_ref*(v(i+1)-v(i))/2 + (v(i+1)+v(i))/2;
end

Z = zeros(numel(V),numel(X));
Z2 = zeros(numel(V),numel(X));
count = 1;
for i=1:num_x
    for j=1:num_v
        coeff = u(count:count+(k+1)^2-1);
        coeff2 = u2(count:count+(k+1)^2-1);
        temp = zeros(poi);
        temp2 = zeros(poi);
        for ii=1:(k+1)
            for jj=1:(k+1)
                temp = temp + coeff((ii-1)*(k+1)+jj)*((leg_vals(jj,:)/sqrt(jac_v))'*(leg_vals(ii,:)/sqrt(jac_x)));
                temp2 = temp2 + coeff2((ii-1)*(k+1)+jj)*((leg_vals(jj,:)/sqrt(jac_v))'*(leg_vals(ii,:)/sqrt(jac_x)));
            end
        end
        Z((j-1)*poi+1:j*poi,(i-1)*poi+1:i*poi) = temp;
        Z2((j-1)*poi+1:j*poi,(i-1)*poi+1:i*poi) = temp2;
        count = count + (k+1)^2;
    end
end
%Z = Z';
Z = flipud(Z);
Z2 = flipud(Z2);
V = flipud(V);

if nargin <= 4
    [XX,VV] = meshgrid(X,V);
    surf(X,V,Z,'EdgeColor','none');
    %caxis([-0.1,1.1]);
    view([0 0 90])
    %scatter3(XX(:),VV(:),Z(:),[],Z(:));
    xlabel(label_x);
    ylabel(label_y);
    title('Plot of numerical soln')
else
    subplot(2,2,1)
    surf(X,V,Z,'LineStyle','none');
    %caxis([-0.1,1.1]);
    view([0 0 90])
    xlabel(label_x);
    ylabel(label_y);
    title('DLRA Solution')
    [XX,VV] = meshgrid(X,V);
    subplot(2,2,2)
    surf(X,V,Z2,'LineStyle','none');
    %caxis([-0.1,1.1]);
    view([0 0 90])
    xlabel(label_x);
    ylabel(label_y);
    title('Full Rank Solution');
    subplot(2,2,3)
    error = Z-Z2;
    surf(XX,VV,error,'EdgeColor','none');
    colorbar
    view([0 0 90])
    xlabel(label_x);
    ylabel(label_y);
    title("Plot of error");
    subplot(2,2,4) 
    if nargin > 5
        plot(myhist(6,:),myhist(2,:));
        title("Rank")
    end
    
end
%surf(X,V,Z-sin(2*pi*X)*(V.^2-1)');
%figure()
%Y = sort(X);
%surf(Y,Y,sin(2*pi*Y)*(Y.^2-1)');

end

