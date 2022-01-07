function [Y] = computeResidualMatrix2(C,S,D,Acell,F,C1,D1)
%Computes the residual (I-P)(LU)*X or ((I-P)(LU))'
%           using low rank decomposition. 
%Here LU = sum A_{i,1}UA_{i,2}' - F.  U = C*S*D' where S is rxr, C,D are nxr
%with r << n.
NN = size(Acell,1);

if nargin < 6
    C1 = C;
    D1 = D;
end

Y = zeros(size(C,1),size(D,1));
% Y = (LU)*X
for l=1:NN
    Y = Y + (Acell{l,1}*C)*S*(Acell{l,2}*D1)';
end
Y = Y - F;
% Projection
%Y = Y - CC'*(LU)*DD'*X
for l=1:NN
    Y = Y - C1*(C1'*Acell{l,1}*C)*S*(D'*Acell{l,2}'*D1)*D1';
end
Y = Y + C1*(C1'*F*D1)*D1';


end


