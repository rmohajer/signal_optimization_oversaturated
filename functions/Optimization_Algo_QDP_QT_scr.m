% This code performs the optimization algorithm for queue discharging
% period (QDP) and using the Queueing theory

%% ARRIVAL FLOWS MAPPING TO FUTURISTIC FLOWS

CHq{1} = ones(N,1)*cmax; % history of cycle times

C = CHq{1};
Time = 0;
Nf = N;
FlowMapping_scr
QHq{1} = qa;

%% CALCULATION OF THE BASIC MATRICES

% the following script calculates the parameters of the delay function for
% optimization.

Delay_params_QT_QD_scr

%% CONSTRAINTS 

% the following script calculates the parameters of the constraints for
% optimization.

Constraints_QT_QD_scr

%% OPTIMIZATION OF THE TOTAL DELAY


LB = THmin;
UB = THmax;
ops = optimoptions('fmincon','Algorithm','sqp'); % active-set OR interior-point-convex

% [x_solN,fval,exitflag] = quadprog(2*AN,BN,H,b,Heq,beq,LB,UB,THmin,ops)

fun_q = @(x) x'*ANq*x+BNq'*x+CNq;
% fun = @(x)BN'*x+CN;
[THETA_q,fval,exitflag] = fmincon(fun_q,THmin,Hq,bq,Heq_q,beq_q,LB,UB,[],ops)



CHq{2} = zeros(N,1);

for j=1:N
    CHq{2}(j) = sum(THETA_q((j-1)*P+1:j*P))+L;
end

C = CHq{2};
Time = 0;
Nf = N;
FlowMapping_scr
QHq{2} = qa;

Delay_params_QT_QD_scr
Constraints_QT_QD_scr

kk=2;
TH0 = THmin;

while (norm(CHq{kk}-CHq{kk-1}))>5/3600*N && exitflag==1 && kk<6

    [THETA_q,fval,exitflag] = fmincon(fun_q,THmax,Hq,bq,Heq_q,beq_q,LB,UB,[],ops);
    
    THETA_q*3600;
    
    disp('end time is')
    sum(C)*3600
    TH0 = THETA_q;
    kk=kk+1;
    CHq{kk} = zeros(N,1);

    for j=1:N
        CHq{kk}(j) = sum(THETA_q((j-1)*P+1:j*P))+L;
    end

    C = CHq{kk};
    sprintf('cycle times:   %d \n',C'*3600)
    Time = 0;
    Nf = N;
    FlowMapping_scr
    QHq{kk} = qa;
    
    disp('qa is: ')
    qa
    
    Delay_params_QT_QD_scr
    Constraints_QT_QD_scr
    norm(QHq{kk}-QHq{kk-1});
    
end

[THETA_q,fval,exitflag] = fmincon(fun_q,THmax,Hq,bq,Heq_q,beq_q,LB,UB,[],ops)
THETA_q*3600
kk = kk+1;
CHq{kk} = zeros(N,1);

for j=1:N
    CHq{kk}(j) = sum(THETA_q((j-1)*P+1:j*P))+L;
end

C = CHq{kk};

% [THETA_q,fval,exitflag] = quadprogIP(2*AN,BN,H,b,Heq,beq,LB,UB)
%  [THETA,fval,exitflag] = cplexqp([],BN,H,b,Heq,beq,LB,UB)
%  THETA*3600

%% TOTAL DELAY AND RESIDUAL QUEUES FOR EACH CYCLE

delta1_q = zeros(P,N);
delta2_q = zeros(P,N);
delta3_q = zeros(P,N);
delta4_q = zeros(P,N);

delay_pk_q = zeros(P,N); % delay of each phase in each cycle
Xj_q = zeros(P,N); % back of the queue in each cycle and each phase

CritTime_q  = zeros(P,3*N);


delta1 = zeros(P,N);
delta2 = zeros(P,N);
delta3 = zeros(P,N);
delta4 = zeros(P,N);

delay_pk = zeros(P,N); % delay of each phase in each cycle
Xj = zeros(P,N); % back of the queue in each cycle and each phase


rownum=1;
for k=1:N

    for p=1:P
        delta1_q(p,1) = delt_0(p)*kj(p);
        delta1(p,1) = delt_0(p);

        % calculation of r_tilde_p,1(k) and r_tilde_p,2(k)
        r_tilpk1 = 0;
        r_tilpk2 = 0;

        if p>1
            for j=1:(p-1)
                r_tilpk1 = r_tilpk1+THETA_q((k-1)*P+j)+l(j);
            end % for j
        end

        if p<P
            for j=(p+1):P
                r_tilpk2 = r_tilpk2+THETA_q((k-1)*P+j)+l(j);
            end
        end
        r_tilpk2 = r_tilpk2+l(p);


        if k>1
            delta1_q(p,k) = delta4_q(p,k-1);
        end

        delta2_q(p,k) = delta1_q(p,k)+qa(p,k)*r_tilpk1;
        delta3_q(p,k) = delta2_q(p,k)-(qc(p)-qa(p,k))*THETA_q((k-1)*P+p);
        delta4_q(p,k) = delta3_q(p,k)+qa(p,k)*r_tilpk2;

        delay_pk_q(p,k) = 0.5*(delta1_q(p,k)+delta2_q(p,k))*r_tilpk1...
            +0.5*(delta3_q(p,k)+delta4_q(p,k))*r_tilpk2;

        Xj_q(p,k) = max([delta2_q(p,k),delta4_q(p,k)]); 
        
        
         if k>1
            delta1(p,k) = delta4(p,k-1);
        end

        delta2(p,k) = delta1(p,k)+Gamma(p,k)*r_tilpk1;
        delta3(p,k) = delta2(p,k)-qc(p)/kj(p)*THETA_q((k-1)*P+p);
        delta4(p,k) = delta3(p,k)+Gamma(p,k)*r_tilpk2;

        delay_pk(p,k) = kj(p)*(delta1(p,k)+0.5*(delta2(p,k)-delta1(p,k)))*r_tilpk1...
            +kj(p)*(delta3(p,k)+0.5*(delta4(p,k)-delta3(p,k)))*r_tilpk2;

        Xj(p,k) = max([delta2(p,k),delta4(p,k)]); 
        
        
        
        rownum = rownum+1;
        
        
        % calculation of the critical points
        if k==1
            CritTime_q(p,(k-1)*4+1) = 0;
            CritTime_q(p,(k-1)*4+2) = CritTime_q(p,(k-1)*4+1)+ r_tilpk1;
            CritTime_q(p,(k-1)*4+3) = CritTime_q(p,(k-1)*4+2)+THETA_q((k-1)*P+p);
            CritTime_q(p,(k-1)*4+4) = CritTime_q(p,(k-1)*4+3)+r_tilpk2;
        else
            CritTime_q(p,(k-1)*4+1) = CritTime_q(p,(k-1)*4);
            CritTime_q(p,(k-1)*4+2) = CritTime_q(p,(k-1)*4+1)+r_tilpk1;
            CritTime_q(p,(k-1)*4+3) = CritTime_q(p,(k-1)*4+2)+THETA_q((k-1)*P+p);
            CritTime_q(p,(k-1)*4+4) = CritTime_q(p,(k-1)*4+3)+r_tilpk2;
        end
        
        


    end % for p
end % for k



TT = sum(C)*60; % total time spent in [min]

if exitflag==1
    % saving the output structure
    info_struct_q(N).delta1_q = delta1_q;
    info_struct_q(N).delta2_q = delta2_q;
    info_struct_q(N).delta3_q = delta3_q;
    info_struct_q(N).delta4_q = delta4_q;
    
    info_struct_q(N).delta1 = delta1;
    info_struct_q(N).delta2 = delta2;
    info_struct_q(N).delta3 = delta3;
    info_struct_q(N).delta4 = delta4;
    
    info_struct_q(N).C = C;
    info_struct_q(N).TT = TT;
    info_struct_q(N).qa = qa;
    info_struct_q(N).delay_q = delay_pk_q;
    info_struct_q(N).Xj_q = Xj_q;
    
    info_struct_q(N).delay = delay_pk;
    info_struct_q(N).Xj = Xj;
    
    info_struct_q(N).Theta = THETA_q;
    info_struct_q(N).exitflag = exitflag;
    info_struct_q(N).CritTime = CritTime_q;
    info_struct_q(N).AN = ANq;
    info_struct_q(N).BN = BNq;

%         figure
%         stairs(1:N,qa(1,:))
%         hold on
%         stairs(1:N,qa(2,:),'r')
%         xlabel('cycle number')
%         ylabel('q^a [veh/h]')
%         legend('Approach 1','Approach 2')
%         fignum = fignum+1;
else
    info_struct_q(N).exitflag = exitflag;
end