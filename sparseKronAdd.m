function [A] = sparseKronAdd(Acell)
%Sparse kron multiply then sum

n = size(Acell,1);
nnzvec = zeros(n,1);
for i=1:n
    nnzvec(i) = nnz(Acell{i,1})*nnz(Acell{i,2});
end
A_nnz = sum(nnzvec);
A_sz = size(Acell{1,2}).*size(Acell{1,1});

I = zeros(A_nnz,1);
J = zeros(A_nnz,1);
S = zeros(A_nnz,1);

count = 1;
for l=1:n
    [i,j,s] = find(kron(Acell{l,2},Acell{l,1}));
    leng = numel(i);
    I(count:count+leng-1) = i;
    J(count:count+leng-1) = j;
    S(count:count+leng-1) = s;
    count = count + leng;
end
I(count:end) = [];
J(count:end) = [];
S(count:end) = [];

A = sparse(I,J,S,A_sz(1),A_sz(2));


end

