function [Acell] = buildSeperableProjectionBlocks(x,v,k,fx,fv)
%Builds the stiffness matrix (Pu,v) where Pu is the projection of u onto
%span{f} where f(x,y)=fx(x)*fy(y).

Fx = buildSeparableSourceX(x,k,fx);
Fv = buildSeparableSourceX(v,k,fv);


%Storing the values of dense matrices, but this can be redone in a
%matrix-free fashion to only store and operator in O(n) instead of O(n^2)
Acell{1,1} = (Fx*Fx')/(Fx'*Fx);
Acell{1,2} = (Fv*Fv')/(Fv'*Fv);
end

