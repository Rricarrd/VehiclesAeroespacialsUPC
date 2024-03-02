function [u,F] = StaticSolver(K,M,TotalDOF,nnodes,g,fixnodes,ndim,Part,Support)

% Generating Neumann and Dirichlett DOFs
[DirichlettDOF, NeumannDOF] = VectorsDOF(TotalDOF, fixnodes, ndim);

%Solves the linear system
%%% K MATRIX SUBPARTS
[KNN, KND, KDN, KDD] = KSubparts(K, NeumannDOF, DirichlettDOF);

%%% FORCES AND DISPLACEMENTS
% Calculation of the Neumann force vector
fN = fNCalc(g,nnodes,M,NeumannDOF); % Force calculation in N

% Calculation of the Dirichlett displacement
if (strcmp(Part,'Part1'))
    uD = zeros(length(DirichlettDOF),1);
elseif(strcmp(Part,'Part2a'))
    uD = zeros(length(DirichlettDOF),1);
    uD(2*(Support-1)*ndim+2)=1e-3; %m to mm %%%%ESTA MALAMENT
elseif(strcmp(Part,'Part2b'))
    
end
% Calculation of the Neumann displacements
uN = KNN\(fN-(KND*uD));

% Calculation of the Dirichlett forces
fD = KDD*uD+KDN*uN;

% Total displacements and forces vectors
u=zeros(TotalDOF,1);
u(DirichlettDOF)=uD;
u(NeumannDOF)=uN;

F=zeros(TotalDOF,1);
F(DirichlettDOF)=fD;
F(NeumannDOF)=fN;

%%% REACTIONS CHECK
% Calculation of the mass in kg
mass = sum(diag(M))/3;
%Check if Dirichlett forces in y direction are equal to weight
w_error = repmat([0 1 0 0 0 0],1,6)*F(DirichlettDOF)+mass*g;
w_error<1e-6;

end