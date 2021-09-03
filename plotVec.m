function [] = plotVec(x,v,k,u,utrue,rel)
%Plots u

if nargin < 6
    rel = false;
end

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

W = zeros(numel(V),numel(X));
Z = zeros(numel(V),numel(X));
count = 1;
for i=1:num_x
    for j=1:num_v
        coeff = u(count:count+(k+1)^2-1);
        temp = zeros(poi);
        for ii=1:(k+1)
            for jj=1:(k+1)
                temp = temp + coeff((ii-1)*(k+1)+jj)*((leg_vals(jj,:)/sqrt(jac_v))'*(leg_vals(ii,:)/sqrt(jac_x)));
            end
        end
        Z((j-1)*poi+1:j*poi,(i-1)*poi+1:i*poi) = temp;
        W((j-1)*poi+1:j*poi,(i-1)*poi+1:i*poi) = w_ref*w_ref';
        count = count + (k+1)^2;
    end
end
%Z = Z';
Z = flipud(Z);
V = flipud(V);

if nargin <= 4
    [XX,VV] = meshgrid(X,V);
    surf(X,V,Z,'EdgeColor','none');
    view([0 0 90])
    %scatter3(XX(:),VV(:),Z(:),[],Z(:));
    xlabel(label_x);
    ylabel(label_y);
    title('Plot of numerical soln')
else
    subplot(2,2,1)
    surf(X,V,Z,'LineStyle','none');
    view([0 0 90])
    xlabel(label_x);
    ylabel(label_y);
    title('Plot of numerical soln')
    [XX,VV] = meshgrid(X,V);
    subplot(2,2,2)
    surf(XX,VV,utrue(XX,VV),'EdgeColor','none');
    view([0 0 90])
    xlabel(label_x);
    ylabel(label_y);
    title('Plot of true soln');
    subplot(2,2,3)
    error = Z-utrue(XX,VV);
    surf(XX,VV,error,'EdgeColor','none');
    colorbar
    view([0 0 90])
    xlabel(label_x);
    ylabel(label_y);
    title("Plot of error");
    subplot(2,2,4)   
    h = findobj(gca,'Type','line');
    T = get(h,'Xdata');
    L2 = get(h,'Ydata');
    if isempty(T)
        T = 1;
        if rel
            L2 = sqrt(sum(error.^2.*W,'all')*jac_x*jac_v)/sqrt(sum(utrue(XX,VV).^2.*W,'all')*jac_x*jac_v);
        else
            L2 = sqrt(sum(error.^2.*W,'all')*jac_x*jac_v);
        end
    else
        T = [T,T(end)+1];
        if rel
            L2 = [L2 sqrt(sum(error.^2.*W,'all')*jac_x*jac_v)/sqrt(sum(utrue(XX,VV).^2.*W,'all')*jac_x*jac_v)];
        else
            L2 = [L2 sqrt(sum(error.^2.*W,'all')*jac_x*jac_v)];
        end
    end
    if numel(T) > 1
        semilogy(T,L2);
    else
        plot(L2,'o');
    end
    if rel
        title("L^{2} Relative Error");
    else
        title("L^{2} Error");
    end
end
%surf(X,V,Z-sin(2*pi*X)*(V.^2-1)');
%figure()
%Y = sort(X);
%surf(Y,Y,sin(2*pi*Y)*(Y.^2-1)');

end

