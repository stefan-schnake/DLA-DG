function [C,S,D] = initHierAdapt(C,S,D,tol,r)
%Using the regularity of the Hierarchcial basis to determine the initial
%number of functions.

N = size(C,1);

norm1C = sum(abs(C));
norm1D = sum(abs(D));

if nargin < 5

    I = find( norm1C < tol*sqrt(N) );
    J = find( norm1D < tol*sqrt(N) );

    fprintf("Hierarchical Adapt: r = %d\n",max([numel(I),numel(J)]));

    if numel(I) < numel(J)
        %Include more vectors to C
        [~,ind] = sort(norm1C);
        I = ind(1:numel(J));
    elseif numel(J) < numel(I)
        %Include more vectors to D
        [~,ind] = sort(norm1D);
        J = ind(1:numel(I));
    end

else
    [~,ind] = sort(norm1C);
    I = ind(1:r);
    [~,ind] = sort(norm1D);
    J = ind(1:r);
end

C = C(:,I);
D = D(:,J);
S = S(I,J);

end

