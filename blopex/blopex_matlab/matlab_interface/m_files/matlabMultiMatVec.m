function blockVectorZ = matlabMultiMatVec(blockVectorX,maskX,...
  blockVectorY,maskY,operator)

%if ~ischar(operator)
if ~isa(operator,'function_handle')
  blockVectorY(:,maskY)=operator*blockVectorX(:,maskX);
else
  blockVectorY(:,maskY) = feval(operator,blockVectorX(:,maskX));
end

blockVectorZ = blockVectorY;