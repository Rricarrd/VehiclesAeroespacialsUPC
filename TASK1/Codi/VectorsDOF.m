function [DirichlettDOF,NeumannDOF] = VectorsDOF(params)
% Generating vector with the fixed DOFs (Dirichlett Conditon)
DirichlettDOF = fixDof(params.fixnodes, params.ndim);

% Genrating vector with the free DOFs (Neumann Conditon)
NeumannDOF = (1:params.totalDOF)';
NeumannDOF(DirichlettDOF) = [];
end

