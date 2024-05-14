% This program investigates total delay minimization for an overflow
% intersection using the dynamic cycle length control method and Shockwave models.

clc
 clear
close('all')
addpath('functions')
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')
addpath('C:\Users\moh044\Dropbox (Sydney Uni)\Transportation\Journal on signal optimization (Oversaturated intersections)\MATLAB codes\quadprogIP-master')

qc2=1440; % veh/h
kj2=140; % veh/km
kc2=65; % veh/km

qc1=2160; % veh/h
kj1=140; % veh/km
kc1=65; % veh/km

qc = [qc1;qc2];
kj = [kj1;kj2];
kc = [kc1;kc2];

%link length
% Del1=1.1; %km
% Del2=1.1;

Del1=.75; %km
Del2=.55;

Del = [Del1;Del2];

cmax = 180/3600; % h

g1min=15/3600; %h
g2min=15/3600;
gmin = [g1min;g2min];


l1=10/3600;l2=10/3600;
l=[l1;l2];
L=l1+l2;

g1max = cmax-L;
g2max = cmax-L;
gmax = [g1max;g2max];

delt1_0=100/kj1; % initial queue length of Phase 1
delt2_0=70/kj2; %initial queue length of Phase 2

% this is for 1 cycle elimination experiment
% delt1_0=60/kj1; 
% delt2_0=40/kj2;


delt_0 = [delt1_0;delt2_0];

T_total = 90/60; %[h]
ts = 3/60;    % sample time for changing flow [h]


N1 = 2; % lower-bound number of cycles that want to discharge
N2 = 10; % upper-bound number of cycles that want to discharge

P = 2; % number of phases
%% GENERATE RANDOM ARRIVAL FLOWS

qmax=[1000, 700];     % veh/h
qmin = [300, 100];% veh/h


% Flow_Generator_scr

load('Arrival_Flows.mat')

%% going to a loop of N
info_struct = struct('delta1',[],'delta2',[],'delta3',[],'delta4',[],'delay',[],...
    'Xj',[],'qa',[],'C',[],'TT',[],'Theta',[],'exitflag',[], 'Time', []);


fignum = 1;
Cycle_numbers_d = [];
flag_delay = 1;
for N=N1:N2
    % the following script runs for optimization procedure
    tstart = tic;
    Optimization_Algo_QDP_LWR_scr
    tend = toc(tstart);
    info_struct(N).solvedur = tend;
    
    if info_struct(N).exitflag>=1
        Cycle_numbers_d = [Cycle_numbers_d, N];
        disp(append(append(int2str(N),', sum is'), int2str((N*L*3600)+sum(info_struct(N).Theta*3600))));
        disp(sum(sum(info_struct(N).delay)))
    end
end


disp('Number of cycles with solution: ')
info_struct.exitflag
disp('The cycle numbers are: ');
disp(Cycle_numbers_d);
disp('The number of feasible cycles are: ');
disp(size(Cycle_numbers_d,2));
info_struct_d = info_struct;

% optimizing the throughput
flag_delay = 0;
Cycle_numbers_th = [];
for N=N1:N2
    % the following script runs for optimization procedure
    tstart = tic;
    Optimization_Algo_QDP_LWR_scr
    tend = toc(tstart);
    info_struct(N).solvedur = tend;
    if info_struct(N).exitflag>=1
        Cycle_numbers_th = [Cycle_numbers_th, N];
        disp(append(append(int2str(N),', sum is'), int2str((N*L*3600)+sum(info_struct(N).Theta*3600))));
        disp(sum(sum(info_struct(N).delay)))
    end
end

disp('Number of cycles with solution: ')
info_struct.exitflag
disp('The cycle numbers are: ');
disp(Cycle_numbers_th);
disp('The number of feasible cycles are: ');
disp(size(Cycle_numbers_th,2));
info_struct_th = info_struct;

%% calculation of the total delay for normal undersaturated operation time.

Nmax = N2;

for flag_delay=[0,1]
    if flag_delay
        info_struct = info_struct_d;
    else
        info_struct = info_struct_th;
    end

    % obtain Nmax in the studied number of cycles to discharge the queues
    if info_struct(Nmax).exitflag==1
        TTmax = info_struct(Nmax).TT;
    else
        while info_struct(Nmax).exitflag ~=1
            Nmax = Nmax-1;
            TTmax = info_struct(Nmax).TT;
        end % while
    end % if

    % plot the arrival flow rate for the lenghtiest scenario
    figure
    stairs(1:Nmax,info_struct(Nmax).qa(1,:))
    hold on
    stairs(1:Nmax,info_struct(Nmax).qa(2,:),'r')
    xlabel('cycle number')
    ylabel('q^a [veh/h]')
    legend('Approach 1','Approach 2')
    fignum = fignum+1;

    % medium cycle time as a threshold to look for additional cycles
    cmed = 1.5; % [min]

    for Np = N1:(Nmax)
        if info_struct(Np).exitflag ==1 && info_struct(Np).TT<TTmax
            % what time is it now?
            Time = sum(info_struct(Np).C); %[h]
            % initial guess on the number of cycles to reach TTmax
            Nt = floor((TTmax-info_struct(Np).TT)/(cmed));
            if Nt>=1
                flag_while = 0;
                while flag_while==0
                    % run the following code for optimization of the signals
                    Optimization_Algo_QDP_LWR_Undersat_scr

                end % while
            else
                info_struct(Np).exitflag_u=0;
            end
        else
            info_struct(Np).exitflag_u=0;
        end
    end

    disp('Number of cycles with solution: ')
    info_struct.exitflag_u
    if flag_delay
        info_struct_d = info_struct;
    else
        info_struct_th = info_struct;
    end
    
end
%% Bang bang algorithm
% [~,dom_approach] = sort(qc,'descend');
% 
% info_structBB = struct('delta1',[],'delta2',[],'delta3',[],'delta4',[],'delay',[],...
%     'Xj',[],'qa',[],'C',[],'TT',[],'Theta',[],'Nopt',[], 'CritTime', []);
% % info_structBB.CritTime = [0,0]';
% % info_structBB.delta3=[.1;.1];
% 
% deltBB_0 = delt_0;
% Time = 0;
% bbloop = 0;
% cmax = 3/60; %[h]
% 
% Res_que0 = diag([0.7,1])*delt_0; % initial residual queue of each phase
% 
% k0 = 0; % this term gets updated to 
% P_larger = 1:P; % the set of phases in which the residual queue size must be larger than zero
% 
% for p_bb=dom_approach'
%     BBexitflag=0;
%     
%     N=1;
%     if bbloop>0
%         if info_structBB.delta3(p_bb,end)<1e-5 
%             BBexitflag=1;
%         end
%     end
%     
%     while BBexitflag==0 && N<=7
%         N=N+1;
%         Optimization_Algo_QDP_LWR_BangBang_scr
%         
%     end
%     P_larger(P_larger==p_bb)=[]; % we remove the constraint that the residual queue of Phase p_bb must be larger than zero
%     for p=1:P
%         deltBB_0(p) = info_structBB.delta4(p,N);
%         Res_que0(p) = info_structBB.delta3(p,N);
%     end
% %     Res_que0(p_bb)=0;
%     
%     Time = sum(info_structBB.C);
%     info_structBB.delta3(p_bb,end);
%     bbloop=bbloop+1;
% end

%% Signal optimization using the queueing theory based appraoch
% 
% info_struct_q = struct('delta1',[],'delta2',[],'delta3',[],'delta4',[],'delay',[],...
%     'Xj',[],'qa',[],'C',[],'TT',[],'Theta',[],'exitflag',[], 'Time', []);
% 
% 
% for N=N1:N2
%     % the following script runs for optimization procedure
%     Optimization_Algo_QDP_QT_scr
% end
% disp('Number of cycles with solution: ')
% info_struct_q.exitflag


%% discharging the queues in one cycle


Optimization_Algo_QDP_LWR_OneCycle_scr
info_struct.exitflag

%% Figures

plotStyle = {'b-','ks-','rx-','g*-','mo-','cd-','yp-'};
plotStyle_com = {'b--','ks--','rx--','g*--','mo--','cd--','yp--'};
% 
Nset = N1:N2;
loopnum = 1;

fignum=1;
% generate figures for comparing the results for different Ns in Nset
Pots_LWR_Ns_scr

%plot the bang bang results and compare with the proposed algorithm

% NBB = sum(info_structBB.Nopt)
% 
% Nset = [NBB-1,NBB,NBB+1];
% loopnum = 1;
% 
% Plots_LWR_BB_scr