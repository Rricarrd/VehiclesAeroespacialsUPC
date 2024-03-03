function [u,F] = StaticSolver(params,mats,nodes,Part,Support)
TotalDOF,nnodes,g,ndim

% Calculation of the Dirichlett displacement
if (strcmp(Part,'Part1'))
    uD = zeros(length(nodes.DirichlettDOF),1);

elseif(strcmp(Part,'Part2a'))
    uD = zeros(length(nodes.DirichlettDOF),1);
    uD((Support-1)*ndim+2)=1e-3; %m to mm 

elseif(strcmp(Part,'Part2b'))
    nodes.DirichlettDOF(2:6:end)=[]; %Remove the y DOF to compute shim height
    nodes.DirichlettDOF(end+1)=(1305-1)*ndim+4; %Prescribed rotation around x
    nodes.DirichlettDOF(end+1)=(1305-1)*ndim+6; %Prescribed rotation around z  
    uD = zeros(length(nodes.DirichlettDOF),1);
    uD(end-1)=500e-6; %Prescribed rotation around x
    uD(end)=-200e-6; %Prescribed rotation around z
    nodes.NeumannDOF = (1:params.totalDOF)';
    nodes.NeumannDOF(nodes.DirichlettDOF) = [];
end

%Solves the linear system
%%% K MATRIX SUBPARTS
[KNN, KND, KDN, KDD] = KSubparts(mats.K, nodes.NeumannDOF, nodes.DirichlettDOF);

%%% FORCES AND DISPLACEMENTS
% Calculation of the Neumann force vector
fN = fNCalc(params.g,params.nnodes,mats.M,nodes.NeumannDOF); % Force calculation in N

% Calculation of the Neumann displacements
uN = KNN\(fN-(KND*uD));

% Calculation of the Dirichlett forces
fD = KDD*uD+KDN*uN;

% Total displacements and forces vectors
u=zeros(params.totalDOF,1);
u(nodes.DirichlettDOF)=uD;
u(nodes.NeumannDOF)=uN;

F=zeros(params.totalDOF,1);
F(nodes.DirichlettDOF)=fD;
F(nodes.NeumannDOF)=fN;



end