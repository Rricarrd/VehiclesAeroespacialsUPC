clear;

%% %%%% PREPROCESSING
%INPUT PARAMETERS
ndim=6;
g=9.81;
fixnodes=[10735 13699 16620 19625 22511 4747];
refnode=1305;
% LOADING DATA
load("fe_model.mat");

% PARAMETERS
TotalDOF=size(K,1);
nnodes = TotalDOF/ndim;

% UNITS
M = M*1e3; % From tons to kg
K = K*1e6;% From kN/mm to N/m




%% %%% CALCULATIONS 
%%% DIRICHLETT AND NEUMANN LOCATIONS
% Generating Neumann and Dirichlett DOFs
[DirichlettDOF, NeumannDOF] = VectorsDOF(TotalDOF, fixnodes, ndim);

%%% K MATRIX SUBPARTS
[KNN, KND, KDN, KDD] = KSubparts(K, NeumannDOF, DirichlettDOF);


%%% FORCES AND DISPLACEMENTS
% Calculation of the Neumann force vector
fN = fNCalc(g, nnodes, M, NeumannDOF); % Force calculation in kN

% Calculation of the Dirichlett displacement
uD = zeros(length(DirichlettDOF),1);

% Calculation of the Neumann displacements

uN = KNN\(fN-(KND*uD));

u=zeros(TotalDOF,1);
u(DirichlettDOF)=uD;
u(NeumannDOF)=uN;

% Calculation of the Dirichlett forces
fD = KDD*uD+KDN*uN;

%%% OTHER CALCULATIONS
% Calculation of the mass in kg
mass = sum(diag(M))/3;


%%%Check if Dirichlett forces in y direction are equal to weight
repmat([0 1 0 0 0 0],1,6)*fD-mass*g<1e-6

%%% Displacements of reference node

u_ref=u((refnode-1)*ndim+1:refnode*ndim)


%%%MODAL ANALYSIS
neig=15;

[MODES EIGENVAL] = eigs(K,M,neig,'sm') ; 
EIGENVAL = diag(EIGENVAL) ; 
FREQ = sqrt(EIGENVAL) ;  

[FREQ,imodes] = sort(FREQ) ;
MODES= MODES(:,imodes); 
