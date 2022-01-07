function [Y] = computeResidual(C,S,D,Acell,F,X,trans,C1,D1)
%Computes the residual (I-P)(LU)*X or ((I-P)(LU))'*X
%           using low rank decomposition. 
%Here LX = sum A_{i,1}XA_{i,2}' - F.  U = C*S*D' where S is rxr, C,D are nxr
%with r << n.
NN = size(Acell,1);

if nargin < 8
    C1 = C;
    D1 = D;
end


if strcmp(trans,'notransp')
    Y = zeros(size(X));
    % Y = (LU)*X
    for l=1:NN
        Y = Y + (Acell{l,1}*C)*S*(Acell{l,2}*D)'*X;
    end
    Y = Y - F*X;
    % Projection is broken into 3 parts
    %Y = Y - LU*(DD')*X
    for l=1:NN
        Y = Y - (Acell{l,1}*C)*S*(D'*Acell{l,2}'*D1)*(D1'*X);
    end
    Y = Y + (F*D1)*(D1'*X);
    %Y = Y + CC'*(LU)*DD'*X
    for l=1:NN
        Y = Y + C1*(C1'*Acell{l,1}*C)*S*(D'*Acell{l,2}'*D1)*(D1'*X);
    end
    Y = Y - C1*(C1'*F*D1)*(D1'*X);
    %Y = Y - CC'*(LU)*X
    for l=1:NN
        Y = Y - C1*(C1'*Acell{l,1}*C)*S*(Acell{l,2}*D)'*X;
    end  
    Y = Y + C1*(C1'*F)*X;
else
    Y = zeros(size(X));
    % Y = (LU)'*X
    for l=1:NN
        Y = Y + (Acell{l,2}*D)*S'*(Acell{l,1}*C)'*X;
    end
    Y = Y - F'*X;
    for l=1:NN
        Y = Y - (Acell{l,2}*D)*S'*(C'*Acell{l,1}'*C1)*(C1'*X);
    end
    Y = Y + (F'*C1)*(C1'*X);
    for l=1:NN
        Y = Y + D1*(D1'*Acell{l,2}*D)*S'*(C'*Acell{l,1}'*C1)*(C1'*X);
    end
    Y = Y - D1*(D1'*F'*C1)*(C1'*X);
    for l=1:NN
        Y = Y - D1*(D1'*Acell{l,2}*D)*S'*(Acell{l,1}*C)'*X;
    end
    Y = Y + D1*(D1'*F')*X;
end
end

