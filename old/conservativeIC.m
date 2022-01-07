function [C,S,D] = conservativeIC(x,v,k,C,S,D,F,G)
%Creates initial conditions that include orthogonal set F/G. 

r = size(S,1);
p = size(F,2);

%Get basis coefficients
f0 = convertMattoVec(x,v,k,C*S*D');

[C,~] = qr([F C],0);
[D,~] = qr([G D],0);

S = zeros(r+p);
for i=1:(r+p)
    for j=1:(r+p)
        S(i,j) = f0'*convertMattoVec(x,v,k,C(:,i)*D(:,j)');
    end
end

end

