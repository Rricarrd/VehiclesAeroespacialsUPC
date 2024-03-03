clear;

%% %%%% PREPROCESSING


% STRUCTURES
params = struct('ndim',[],'g',[],'fixnodes',[],'refnode',[],'totalDOF',[],'nnodes',[]);
nodes = struct('DirichlettDOF',[],'NeumannDOF',[]);
mats = struct('K',[],'M',[]);
submats = struct('KNN',[],'KND',[],'KDN',[],'KDD',[],'MNN',[]);
results = struct('u',[],'f',[],'w_error',[],'mass',[]);

%INPUT PARAMETERS
params.ndim=6;
params.g=9.81;
params.fixnodes=[10735 13699 16620 19625 22511 4747];
params.refnode=1305;

% LOADING DATA
load("fe_model.mat");

% PARAMETERS
params.totalDOF=size(K,1);
params.nnodes = params.totalDOF/params.ndim;

% UNITS

mats.M = M*1000; % From tons to kg
mats.K = K*1000; % From N/mm to N/m



%%%%%% STATIC ANALYSIS
%% %%% PART 1
%%% DIRICHLETT AND NEUMANN LOCATIONS

% Generating Neumann and Dirichlett DOFs
[nodes.DirichlettDOF, nodes.NeumannDOF] = VectorsDOF(params);

% BODY FORCES SOLUTION
Part='Part1';
Support=1;
[results.u,results.F] = StaticSolver(params,mats,nodes,Part,Support);

%%% REACTIONS CHECK
% Calculation of the mass in kg
results.mass = sum(diag(M))/3;
%Check if Dirichlett forces in y direction are equal to weight
results.w_error = repmat([0 1 0 0 0 0],1,params.ndim)*results.F(nodes.DirichlettDOF)+results.mass*params.g;

%PRINT RESULTS TO HDF5
uhdf=zeros(nnodes,ndim);
for i=1:ndim
    uhdf(:,i)=u(i:ndim:end);
end
fillhdf("template.h5","Part1Output.h5",uhdf);

%% %%%% PART 2
%%% DISPLACEMENTS AND ROTATIONS
%Compute the displacements and rotations of the referencenode due
%to enforced unit displacements in y direction in each support.

% Displacements of reference node
u_ref=u((refnode-1)*ndim+1:refnode*ndim);

%Change of prescribed displacements 
Support=1;
Part='Part2a';
[u,F] = StaticSolver(K,M,totalDOF,nnodes,g,DirichlettDOF,NeumannDOF,ndim,Part,Support);


% SHIMS CALCULATIONS
fixnodes=[10735 13699 16620 19625 22511 4747 1305];
Part='Part2b';
[u_shim,F_shim] = StaticSolver(K,M,totalDOF,nnodes,g,DirichlettDOF,NeumannDOF,ndim,Part,Support);

u_shim=u_shim(DirichlettDOF);
u_shim=u_shim(2:6:end);

%%%%%% DYNAMIC ANALYSIS
%% %%% PART 1: MODAL ANALYSIS
neig=11;
f_max = 2000;
[KNN,KND, KDN, KDD] = KSubparts(K, NeumannDOF, DirichlettDOF);

% Constrained masses
MNN = M(NeumannDOF,NeumannDOF);

% Modal response constrained
[MODES, FREQ] = modResponse(KNN,MNN,neig);

% Modal response unconstrained
[MODES_U, FREQ_U] = modResponse(K,M,neig);

%% %%% PART 2: RESPONSE ANALYSIS
damping = 0;
values = 500;
% Unconstrained nodes excitation in X direction at 1m/s^2
[X, ~] = freqResponse(values,f_max,nnodes,NeumannDOF,MODES,M,MNN,KNN,damping,FREQ);
X_ref = X(refnode,:);

% Unconstrained nodes excitation in X direction at 1m/s^2
damping = 0.02;
[X_damped, f_values] = freqResponse(values,f_max,nnodes,NeumannDOF,MODES,M,MNN,KNN,damping,FREQ);
X_ref_damped = X_damped(refnode,:);

figure(1)
plot(f_values, abs(X_ref))

hold on
plot(f_values, abs(X_ref_damped))
hold off
