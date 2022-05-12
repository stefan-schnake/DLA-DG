function [U] = steadyStateSoln(x,v,k,Acell)
%Solves steady state problem Acell

A = zeros(numel(Acell{1,1}));
for i=1:size(Acell,1)
    A = A + kron(Acell{i,2},Acell{i,1});
end


end

