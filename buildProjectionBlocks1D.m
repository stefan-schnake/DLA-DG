function [Acell] = buildProjectionBlocks1D(x,y,k,fy,alpha)
%Builds the stiffness matrix alpha(u-f_yPu,v) where Pu is the projection of u onto
%span{f_y} w.r.t \W_y

Fy = buildSeparableSourceX(y,k,fy);


%Storing the values of dense matrices, but this can be redone in a
%matrix-free fashion to only store and operator in O(n) instead of O(n^2)
Acell{1,1} = speye((numel(x)-1)*(k+1));
Acell{1,2} = alpha*(speye((numel(y)-1)*(k+1))-(Fy*Fy')/(Fy'*Fy));
end