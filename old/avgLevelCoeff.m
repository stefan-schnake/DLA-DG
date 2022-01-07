function [avg] = avgLevelCoeff(k,coeffs)
%Determines the average abs coefficient value on each sublevel in a
%hierarchical basis
[n,r] = size(coeffs);

m = log2(n/(k+1))+1;

ind = zeros(m+1,1);
ind(1) = 1;
ind(2) = k+2;
for j=2:m
    ind(j+1) = ind(j) + (k+1)*2^(j-2);
end

avg = zeros(m,r);
%Populate first row
for i=1:r
    for j=1:m
        avg(j,i) = norm(coeffs(ind(j):(ind(j+1)-1),i));
    end
end

end

