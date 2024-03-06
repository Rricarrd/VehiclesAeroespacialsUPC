clear;
clc;
%% %%%% PREPROCESSING
%INPUT PARAMETERS
ndim=6;
nsuports=6;
g=9.81;
fixnodes=[10735 13699 16620 19625 22511 4747];
refnode=1305;
LATEX=0; %Parameter that controls the output for inserting in the LaTeX document
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

if(LATEX)
    Faux=F(DirichlettDOF);
    Fresults=zeros(nsuports,ndim);
    for i=1:ndim
        Fresults(:,i)=Faux(i:ndim:end); 
    end
    for j=1:nsuports
        for i=1:ndim-1
            fprintf('%.2f & \n',Fresults(j,i))
        end
        fprintf('%.2f \\\\ \\hline \n',Fresults(j,i+1))
        fprintf('\\textbf{Suport %d} & \n',j+1)
    end
    Fcheck(1,:)=sum(Fresults,1);
end 
%%% REACTIONS CHECK
[w_error,mass] = reactionsCheck(M,F,g,DirichlettDOF,dir);


% Displacements of reference node
u_ref=u((refnode-1)*ndim+1:refnode*ndim)*10^6  % back to um;

% Forces at supports
F_supports = F(DirichlettDOF);

%PRINT RESULTS TO HDF5
uhdf=zeros(nnodes,ndim);
for i=1:ndim
    uhdf(:,i)=u(i:ndim:end)*1000; % back to mm
end
output_name = "Output/Part1Output" + dir + ".h5";
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

% Displacements of reference node after imposing displacements of the
% dirichlett nodes in y direction
u_ref_2a=u((refnode-1)*ndim+1:refnode*ndim)*1e3; % m to mm

% SHIMS CALCULATIONS
fixnodes=[10735 13699 16620 19625 22511 4747 1305];
Part='Part2b';
[u_shim,F_shim] = StaticSolver(K,M,TotalDOF,nnodes,g,DirichlettDOF,NeumannDOF,ndim,Part,Support,dir);

u_shim=u_shim(DirichlettDOF);
u_shim=u_shim(2:6:end);
shim_size=1+u_shim;

Part='ShimCheck';
[u_check,F_check] = StaticSolver(K,M,TotalDOF,nnodes,g,DirichlettDOF,NeumannDOF,ndim,Part,Support,dir,shim_size);
u_ref_check=u_check((refnode-1)*ndim+1:refnode*ndim)*1e+6;

%%%%%% DYNAMIC ANALYSIS
%% %%% PART 1: MODAL ANALYSIS
neig=11;
f_max = 2000;
[KNN,KND, KDN, KDD] = KSubparts(K, NeumannDOF, DirichlettDOF);

% Constrained masses
MNN = M(NeumannDOF,NeumannDOF);

%% Modal response constrained
[MODES, FREQ] = modResponse(KNN,MNN,neig);
FREQ_HZ = FREQ/(2*pi);

% Build the displacements for each node
MODES_TOT=zeros(TotalDOF,neig);
for i = 1:neig
    MODES_TOT(DirichlettDOF, i) = 0;
    MODES_TOT(NeumannDOF, i)=MODES(:,i);
end

% Save modes in the uhdf file
for j = 1:neig-6
    uhdf=zeros(nnodes,ndim);
    for i=1:ndim
        column = MODES_TOT(:, j);
        uhdf(:,i)=column(i:ndim:end)*1000; % back to mm
    end
    output_name = "Output/Part3OutputConstrained" + j + ".h5";
    fillhdf("template.h5",output_name,uhdf);
end

%% Modal response unconstrained
[MODES_U, FREQ_U] = modResponse(K,M,neig);
FREQ_HZ_U = FREQ_U/(2*pi);

% Build the displacements for each node
MODES_TOT_U=zeros(TotalDOF,neig);
for i = 1:neig
    MODES_TOT_U(1:TotalDOF, i)=MODES_U(:,i);
end

% Save modes in the uhdf file
for j = 7:neig
    uhdf=zeros(nnodes,ndim);
    for i=1:ndim
        column = MODES_TOT_U(:, j);
        uhdf(:,i)=column(i:ndim:end)*1000; % back to mm
    end
    output_name = "Output/Part3OutputUnconstrained" + j + ".h5";
    fillhdf("template.h5",output_name,uhdf);
end


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
% plot(f_values(1:values/2), abs(X_ref(1:values/2))*10^6)

hold on
plot(f_values(1:values/2), abs(X_ref_damped(1:values/2))*10^6)
hold off

xlabel("Frequency (Hz)", Interpreter="latex")
ylabel("Amplitude ($\mu m$)", Interpreter="latex")
title("Damped frequency response up to 2000Hz", Interpreter="latex")
