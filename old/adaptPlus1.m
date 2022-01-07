function [C,S,D] = adaptPlus1(C,S,D)
%Adds a column to C and S that is maintaing orthogonality
%Adds 1 zero row and column to S

q = createNewVec(C);
while norm(q) < 1e-13
    fprintf('Vector chosen was in range of matrix, retrying...\n');
end
q = q/norm(q);
C = [C q];

q = createNewVec(D);
while norm(q) < 1e-13
    fprintf('Vector chosen was in range of matrix, retrying...\n');
end
q = q/norm(q);
D = [D q];

S = [S zeros(size(S,1),1)];
S = [S;zeros(size(S,2),1)'];

end

function q = createNewVec(C)
    q = rand(size(C,1),1);
    q = q - C*C'*q;
end

