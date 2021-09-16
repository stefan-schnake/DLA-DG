function [FU] = computeFullRank(U,Acell,BC)
%Computes full rank F(U)

FU = zeros(size(U));
for i=1:size(Acell,1)
    FU = FU + Acell{i,1}*U*Acell{i,2}';
end
if BC.use
    for l=1:size(BC.cell1,1)
        FU = FU + BC.cell1{l,1}*BC.cell1{l,2}*BC.cell1{l,3}';
    end
end
end

