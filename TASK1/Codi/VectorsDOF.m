function [DirichlettDOF,NeumannDOF] = VectorsDOF(TotalDOF, fixnodes, ndim)
% Generating vector with the fixed DOFs (Dirichlett Conditon)
DirichlettDOF = fixDof(fixnodes, ndim);

% Genrating vector with the free DOFs (Neumann Conditon)
NeumannDOF = (1:TotalDOF)';
NeumannDOF(DirichlettDOF) = [];
end

