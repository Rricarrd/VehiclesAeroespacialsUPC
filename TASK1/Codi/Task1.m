clear;

%% %%%% PREPROCESSING
%INPUT PARAMETERS
ndim=6;
g=9.81;
fixnodes=[10735 13699 16620 19625 22511 4747];

% LOADING DATA
load("fe_model.mat");

% PARAMETERS
TotalDOF=size(K,1);
nnodes = TotalDOF/ndim;

% UNITS
M = M*1e3; % From tons to kgs
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

% Calculation of the Dirichlett forces
fD = KDD*uD+KDN*uN;

%%% OTHER CALCULATIONS
% Calculation of the mass in kg
W = sum(fN)/-9.81;

