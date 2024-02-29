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
K = K*1e3; % From N/mm to N/m


%% %%% STATIC ANALYSIS
%%% DIRICHLETT AND NEUMANN LOCATIONS
% Generating Neumann and Dirichlett DOFs
[DirichlettDOF, NeumannDOF] = VectorsDOF(TotalDOF, fixnodes, ndim);

[u,F] = StaticSolver(K,M,NeumannDOF,DirichlettDOF,nnodes);


%%% REACTIONS CHECK
% Calculation of the mass in kg
mass = sum(diag(M))/3;
%Check if Dirichlett forces in y direction are equal to weight
w_error = repmat([0 1 0 0 0 0],1,6)*fD+mass*g;
w_error<1e-6;
%%% Displacements of reference node
u_ref=u((refnode-1)*ndim+1:refnode*ndim);

%PRINT RESULTS TO HDF5
uhdf=zeros(nnodes,ndim);
for i=1:ndim
    uhdf(:,i)=u(i:ndim:end);
end
fillhdf("template.h5","Part1Output.h5",uhdf);

%% %%%%
%PART 2
%Compute the displacements and rotations of 
% the reference node due to enforced unit displacements in y direction in
%each support.
%Change of prescribed displacements 
uD = zeros(length(DirichlettDOF),1);
Support=1;
uD(2*(Support-1)*ndim+1)=1e-3; %m to mm

[u,F] = StaticSolver(K,M,NeumannDOF,DirichlettDOF,nnodes);



%% %%% DYNAMIC ANALYSIS

%%% MODAL ANALYSIS
neig=10;

% Constrained masses
MNN = M(NeumannDOF,NeumannDOF);

% Modal response constrained
[MODES, FREQ] = modResponse(KNN,MNN,neig);

% Modal response unconstrained
[MODES_U, FREQ_U] = modResponse(K,M,neig);

%%% RESPONSE ANALYSIS



