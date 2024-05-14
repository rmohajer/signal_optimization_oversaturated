% This program investigates total delay minimization for an overflow
% intersection using the dynamic cycle length control method and Shockwave models.

clc
 clear
close('all')
addpath('functions')
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')
addpath('C:\Users\moh044\Dropbox (Sydney Uni)\Transportation\Journal on signal optimization (Oversaturated intersections)\MATLAB codes\quadprogIP-master')

qc2=1500;
kj2=140;
kc2=65;

qc1=1500;
kj1=140;
kc1=65;


qc = [qc1;qc2];
kj = [kj1;kj2];
kc = [kc1;kc2];

%link length
Del1=.9; %km
Del2=.9;

% Del1=1.0; %km
% Del2=1.0;

Del = [Del1;Del2];

cmax = 400/3600; % h

g1min=15/3600; %h
g2min=15/3600;
gmin = [g1min;g2min];


l1=7/3600;l2=7/3600;
l=[l1;l2];
L=l1+l2;

g1max = cmax-L;
g2max = cmax-L;
gmax = [g1max;g2max];


delt1_0=8/kj1; 
delt2_0=0/kj2;

delt_0 = [delt1_0;delt2_0];

%Congestion period:
T_total = 40/60; %[h]
ts = 10/60;    % sample time for changing flow [h]


N1 = floor(T_total/cmax); % lower-bound number of cycles that want to discharge
N2 = 18; % upper-bound number of cycles that want to discharge

P = 2; % number of phases

% alphas for each phase
global Alphas
Alphas = [1.3,1.3];

omega1=.3;
omega2=.7;
%% GENERATE RANDOM ARRIVAL FLOWS

qmax=[800, 1000];     % veh/h
qmin = [400, 400];% veh/h
T_generate = 60/60; %[h]

% Flow_Generator_QF_scr

load('Arrival_Flows_QFP.mat')

% Qa = [700*ones(1,5);900*ones(1,5)];


%% going to a loop of N and optimization using the total delay function as objective
info_struct_QF = struct('delta1',[],'delta2',[],'delta3',[],'delta4',[],'delay',[],...
    'Xj',[],'qa',[],'C',[],'TT',[],'Theta',[],'exitflag',[], 'Time', []);

Cycle_numbers_d = [];
flag_delay = 1;
for N=N1:N2
    % the following script runs for optimization procedure
    Optimization_Algo_QFP_scr
    if info_struct_QF(N).exitflag>=1
        Cycle_numbers_d = [Cycle_numbers_d, N];
        disp(append(append(int2str(N),', sum is'), int2str((N*L*3600)+sum(info_struct_QF(N).Theta*3600))));
        disp(sum(sum(info_struct_QF(N).delay)))
    end
end

disp('Number of cycles with solution: ')
info_struct_QF.exitflag
disp('The cycle numbers are: ');
disp(Cycle_numbers_d);
disp('The number of feasible cycles are: ');
disp(size(Cycle_numbers_d,2));
info_struct_QF_d = info_struct_QF;

% optimizing the throughput
flag_delay = 2;
Cycle_numbers_th = [];
for N=N1:N2
    % the following script runs for optimization procedure
    Optimization_Algo_QFP_scr
    
    if info_struct_QF(N).exitflag>=1
        Cycle_numbers_th = [Cycle_numbers_th, N];
        disp(append(append(int2str(N),', sum is'), int2str((N*L*3600)+sum(info_struct_QF(N).Theta*3600))));
        disp(sum(sum(info_struct_QF(N).delay)))
    end
end

disp('Number of cycles with solution: ')
info_struct_QF.exitflag
disp('The cycle numbers are: ');
disp(Cycle_numbers_th);
disp('The number of feasible cycles are: ');
disp(size(Cycle_numbers_th,2));
info_struct_QF_th = info_struct_QF;



% optimizing the mixed objective
flag_delay = 3;
Cycle_numbers_sp = [];
for N=N1:N2
    % the following script runs for optimization procedure
    Optimization_Algo_QFP_scr
    
    if info_struct_QF(N).exitflag>=1
        Cycle_numbers_sp = [Cycle_numbers_sp, N];
        disp(append(append(int2str(N),', sum is'), int2str((N*L*3600)+sum(info_struct_QF(N).Theta*3600))));
        disp(sum(sum(info_struct_QF(N).delay)))
    end
end

disp('Number of cycles with solution: ')
info_struct_QF.exitflag
disp('The cycle numbers are: ');
disp(Cycle_numbers_sp);
disp('The number of feasible cycles are: ');
disp(size(Cycle_numbers_sp,2));
info_struct_QF_sp = info_struct_QF;



%% Figures

plotStyle = {'b-','ks-','rx-','g*-','mo-','cd-','yp-','k--'};
plotStyle_com = {'b--','ks--','rx--','g*--','mo--','cd--','yp--','k--'};
% 
% Nset = [10,12, 14,16,18,20];
Nset = 11:N2;

loopnum = 1;

fignum=1;
% generate figures for comparing the results for different Ns in Nset
Pots_LWR_QFP_scr