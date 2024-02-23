function fixDOF = fixDof(fixnodes,ndim)
% Fucntion that generates a list with the fixed nodes od the 
%
% Args:
%
% Returns:
%

fixDOF=zeros(length(fixnodes),ndim);
for i=1:length(fixnodes)
    fixDOF(i,:) =(fixnodes(i)-1)*ndim+1:fixnodes(i)*ndim;
end
fixDOF=reshape(fixDOF',1,[])';

end

