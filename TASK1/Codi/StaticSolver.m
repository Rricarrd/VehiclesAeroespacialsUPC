function [u,F] = StaticSolver(K,M,NeumannDOF,DirichlettDOF,TotalDOF,nnodes,g)
%Solves the linear system
%%% K MATRIX SUBPARTS
[KNN, KND, KDN, KDD] = KSubparts(K, NeumannDOF, DirichlettDOF);

%%% FORCES AND DISPLACEMENTS
% Calculation of the Neumann force vector
fN = fNCalc(g,nnodes,M,NeumannDOF); % Force calculation in N

% Calculation of the Dirichlett displacement
uD = zeros(length(DirichlettDOF),1);

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


end