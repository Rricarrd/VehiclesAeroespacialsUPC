function [u,F] = StaticSolver(K,M,TotalDOF,nnodes,g,DirichlettDOF,NeumannDOF,ndim,Part,Support,dir,shim_size)

% Calculation of the Dirichlett displacement
if (strcmp(Part,'Part1'))
    uD = zeros(length(DirichlettDOF),1);

elseif(strcmp(Part,'Part2a'))
    uD = zeros(length(DirichlettDOF),1);
    uD((Support-1)*ndim+2)=1e-3; %m to mm 

elseif(strcmp(Part,'Part2b'))
    % DirichlettDOF(2:6:end)=[]; %Remove the y DOF to compute shim height
    % DirichlettDOF(end+1)=(1305-1)*ndim+4; %Prescribed rotation around x
    % DirichlettDOF(end+1)=(1305-1)*ndim+6; %Prescribed rotation around z  
    % uD = zeros(length(DirichlettDOF),1);
    % uD(end-1)=500e-6; %Prescribed rotation around x
    % uD(end)=-200e-6; %Prescribed rotation around z
    % NeumannDOF = (1:TotalDOF)';
    % NeumannDOF(DirichlettDOF) = [];
    DirichlettDOF(2:6:end)=[]; %Remove the y DOF to compute shim height
    DirichlettDOF=vertcat(DirichlettDOF,(1305-1)*ndim+[1:6]'); %Prescribed rotation around x
    uD = zeros(length(DirichlettDOF),1);
    uD(end-2)=500e-6; %Prescribed rotation around x
    uD(end)=-200e-6; %Prescribed rotation around z
    NeumannDOF = (1:TotalDOF)';
    NeumannDOF(DirichlettDOF) = [];
elseif(strcmp(Part,'ShimCheck'))
    uD = zeros(length(DirichlettDOF),1);
      %uD(2:6:end)=1e-3+[0.0333 0.0602 0.0289 -0.0432 -0.0621 -0.0172]; %ushim-mean(ushim)
      uD(2:6:end)=shim_size;
        %1e-3+[0.0074 0.0432 0.0438 -0.0515 -0.0512 0.0083]

end

%Solves the linear system
%%% K MATRIX SUBPARTS
[KNN, KND, KDN, KDD] = KSubparts(K, NeumannDOF, DirichlettDOF);

%%% FORCES AND DISPLACEMENTS
% Calculation of the Neumann force vector
fN = fNCalc(g,nnodes,M,NeumannDOF, dir); % Force calculation in N

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