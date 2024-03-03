clear;
clc;
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

M = M*1000; % From tons to kg
K = K*1000; % From N/mm to N/m


%%%%%% STATIC ANALYSIS
%% %%% PART 1
%%% DIRICHLETT AND NEUMANN LOCATIONS

% Generating Neumann and Dirichlett DOFs
[DirichlettDOF, NeumannDOF] = VectorsDOF(TotalDOF, fixnodes, ndim);

% BODY FORCES SOLUTION
Part='Part1';
Support=1;
dir = 'Z';
[u,F] = StaticSolver(K,M,TotalDOF,nnodes,g,DirichlettDOF,NeumannDOF,ndim,Part,Support,dir);

%%% REACTIONS CHECK
[w_error,mass] = reactionsCheck(M,F,g,DirichlettDOF,dir);

% Displacements of reference node
u_ref=u((refnode-1)*ndim+1:refnode*ndim)*10^6;

% Forces at supports
F_supports = F(DirichlettDOF);

%PRINT RESULTS TO HDF5
uhdf=zeros(nnodes,ndim);
for i=1:ndim
    uhdf(:,i)=u(i:ndim:end);
end
output_name = "Part1Output" + dir + ".h5";
fillhdf("template.h5",output_name,uhdf);

%% %%%% PART 2
%%% DISPLACEMENTS AND ROTATIONS
%Compute the displacements and rotations of the referencenode due
%to enforced unit displacements in y direction in each support.


%Change of prescribed displacements 
dir = 'Y';
Support=1;
Part='Part2a';
[u,F] = StaticSolver(K,M,TotalDOF,nnodes,g,DirichlettDOF,NeumannDOF,ndim,Part,Support,dir);


% SHIMS CALCULATIONS
fixnodes=[10735 13699 16620 19625 22511 4747 1305];
Part='Part2b';
[u_shim,F_shim] = StaticSolver(K,M,TotalDOF,nnodes,g,DirichlettDOF,NeumannDOF,ndim,Part,Support,dir);

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
