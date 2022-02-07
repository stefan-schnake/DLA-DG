function [] = plotVecX(x,k,u)
%Plots u

num_x = numel(x)-1;
poi = 7;

[quad_ref, ~]  = lgwt(poi,-1,1);
quad_ref = quad_ref';
quad_ref = fliplr(quad_ref);

[leg_vals,~,~,~] = buildLegendre(poi,k);
leg_vals = fliplr(leg_vals);
leg_vals = leg_vals/sqrt((x(2)-x(1))/2);

X = zeros(poi*num_x,1);
for i=1:num_x
    X((i-1)*poi+1:i*poi) = quad_ref*(x(i+1)-x(i))/2 + (x(i+1)+x(i))/2;
end

Z = zeros(size(X));
count = 1;
for i=1:num_x
        coeff = u(count:count+(k+1)-1);
        temp = leg_vals'*coeff;
        Z((i-1)*poi+1:i*poi) = temp;
        count = count + (k+1);
end

plot(X,Z);
%surf(X,V,Z-sin(2*pi*X)*(V.^2-1)');
%figure()
%Y = sort(X);
%surf(Y,Y,sin(2*pi*Y)*(Y.^2-1)');

end

