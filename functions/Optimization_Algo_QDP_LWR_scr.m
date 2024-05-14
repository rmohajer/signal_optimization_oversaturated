% This code performs the optimization algorithm for queue discharging
% period and using the LWR theory

%% ARRIVAL FLOWS MAPPING TO FUTURISTIC FLOWS

CH{1} = ones(N,1)*cmax; % history of cycle times

C = CH{1};
Time = 0;
Nf = N;
FlowMapping_scr
QH{1} = qa;

%% CALCULATION OF THE BASIC MATRICES

% the following script calculates the parameters of the delay function for
% optimization.

Delay_params_LWR_QD_scr


%% CONSTRAINTS 

% the following script calculates the parameters of the constraints for
% optimization.

Constraints_LWR_QD_scr

%% OPTIMIZATION OF THE TOTAL DELAY


LB = THmin;
UB = THmax;
ops = optimoptions('fmincon','Algorithm','sqp'); % active-set OR interior-point-convex

% [x_solN,fval,exitflag] = quadprog(2*AN,BN,H,b,Heq,beq,LB,UB,THmin,ops)

if flag_delay
    fun = @(x) x'*AN*x+BN'*x+CN;
else
    fun = @(x) -repmat(qc,N,1)'*x;  % for throughput minimization
end

[THETA,fval,exitflag] = fmincon(fun,THmax,H,b,Heq,beq,LB,UB,[],ops)



CH{2} = zeros(N,1);

for j=1:N
    CH{2}(j) = sum(THETA((j-1)*P+1:j*P))+L;
end

C = CH{2};
Time = 0;
Nf = N;
FlowMapping_scr
QH{2} = qa;

Delay_params_LWR_QD_scr
Constraints_LWR_QD_scr

kk=2;
TH0 = THmin;

while (norm(QH{kk}-QH{kk-1}))> 200 && exitflag>=1 && kk<6

    [THETA,fval,exitflag] = fmincon(fun,TH0,H,b,Heq,beq,LB,UB,[],ops);
    
    THETA*3600;
    
    disp('end time is')
    sum(C)*3600
    TH0 = THETA;
    kk=kk+1;
    CH{kk} = zeros(N,1);

    for j=1:N
        CH{kk}(j) = sum(THETA((j-1)*P+1:j*P))+L;
    end

    C = CH{kk};
    sprintf('cycle times:   %d \n',C'*3600)
    Time = 0;
    Nf = N;
    FlowMapping_scr
    QH{kk} = qa;
    
    disp('qa is: ')
    qa
    
    Delay_params_LWR_QD_scr
    Constraints_LWR_QD_scr
    norm(QH{kk}-QH{kk-1});
    
end

[THETA,fval,exitflag] = fmincon(fun,TH0,H,b,Heq,beq,LB,UB,[],ops)
THETA*3600
kk = kk+1;
CH{kk} = zeros(N,1);

for j=1:N
    CH{kk}(j) = sum(THETA((j-1)*P+1:j*P))+L;
end

C = CH{kk};

% [THETA,fval,exitflag] = quadprogIP(2*AN,BN,H,b,Heq,beq,LB,UB)
%  [THETA,fval,exitflag] = cplexqp([],BN,H,b,Heq,beq,LB,UB)
%  THETA*3600
%% TOTAL DELAY AND RESIDUAL QUEUES FOR EACH CYCLE

delta1 = zeros(P,N);
delta2 = zeros(P,N);
delta3 = zeros(P,N);
delta4 = zeros(P,N);

delay_pk = zeros(P,N); % delay of each phase in each cycle
TH_pk = zeros(P,N); % throughput of each phase in each cycle
Xj = zeros(P,N); % back of the queue in each cycle and each phase

CritTime  = zeros(P,3*N);


rownum=1;
for k=1:N

    for p=1:P
        delta1(p,1) = delt_0(p);

        % calculation of r_tilde_p,1(k) and r_tilde_p,2(k)
        r_tilpk1 = 0;
        r_tilpk2 = 0;

        if p>1
            for j=1:(p-1)
                r_tilpk1 = r_tilpk1+THETA((k-1)*P+j)+l(j);
            end % for j
        end

        if p<P
            for j=(p+1):P
                r_tilpk2 = r_tilpk2+THETA((k-1)*P+j)+l(j);
            end
        end
        r_tilpk2 = r_tilpk2+l(p);


        if k>1
            delta1(p,k) = delta4(p,k-1);
        end

        delta2(p,k) = delta1(p,k)+Gamma(p,k)*r_tilpk1;
        delta3(p,k) = delta2(p,k)-qc(p)/kj(p)*THETA((k-1)*P+p);
        delta4(p,k) = delta3(p,k)+Gamma(p,k)*r_tilpk2;

        delay_pk(p,k) = kj(p)*(delta1(p,k)+0.5*(delta2(p,k)-delta1(p,k)))*r_tilpk1...
            +kj(p)*(delta3(p,k)+0.5*(delta4(p,k)-delta3(p,k)))*r_tilpk2;
        
        TH_pk(p,k) = qc(p)*THETA((k-1)*P+p);

        Xj(p,k) = max([delta2(p,k),delta4(p,k)]); 
        rownum = rownum+1;

        % calculation of the critical points
        if k==1
            CritTime(p,(k-1)*4+1) = 0;
            CritTime(p,(k-1)*4+2) = CritTime(p,(k-1)*4+1)+ r_tilpk1;
            CritTime(p,(k-1)*4+3) = CritTime(p,(k-1)*4+2)+THETA((k-1)*P+p);
            CritTime(p,(k-1)*4+4) = CritTime(p,(k-1)*4+3)+r_tilpk2;
        else
            CritTime(p,(k-1)*4+1) = CritTime(p,(k-1)*4);
            CritTime(p,(k-1)*4+2) = CritTime(p,(k-1)*4+1)+r_tilpk1;
            CritTime(p,(k-1)*4+3) = CritTime(p,(k-1)*4+2)+THETA((k-1)*P+p);
            CritTime(p,(k-1)*4+4) = CritTime(p,(k-1)*4+3)+r_tilpk2;
        end
        
        


    end % for p
end % for k



TT = sum(C)*60; % total time spent in [min]

if exitflag>=1
    % saving the output structure
    info_struct(N).delta1 = delta1;
    info_struct(N).delta2 = delta2;
    info_struct(N).delta3 = delta3;
    info_struct(N).delta4 = delta4;
    info_struct(N).C = C;
    info_struct(N).TT = TT;
    info_struct(N).qa = qa;
    info_struct(N).delay = delay_pk;
    info_struct(N).through = TH_pk;
    info_struct(N).Xj = Xj;
    info_struct(N).Theta = THETA;
    info_struct(N).exitflag = exitflag;
    info_struct(N).CritTime = CritTime;
    info_struct(N).AN = AN;
    info_struct(N).BN = BN;
    info_struct(N).Eta = Eta;

else
    info_struct(N).exitflag = exitflag;
end