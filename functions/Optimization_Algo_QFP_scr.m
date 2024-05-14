% This code performs the optimization algorithm for queue discharging
% period and using the LWR theory

%% ARRIVAL FLOWS MAPPING TO FUTURISTIC FLOWS

CH{1} = ones(N,1)*T_total/N; % history of cycle times

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

Constraints_LWR_QFP_scr

%% OPTIMIZATION OF THE TOTAL DELAY


LB = THmin;
UB = THmax;
ops = optimoptions('fmincon','Algorithm','sqp'); % active-set OR interior-point

% [x_solN,fval,exitflag] = quadprog(2*AN,BN,H,b,Heq,beq,LB,UB,THmin,ops)

if flag_delay==1
    fun = @(x) x'*AN*x+BN'*x+CN;
elseif flag_delay==2
    fun = @(x) -repmat(qc,N,1)'*x;  % for throughput minimization
else
    fun = @(x) omega1*(x'*AN*x+BN'*x+CN)+0.5*omega2*SP_fun(x,Del,delt_0,Gamma,qc,kj,l,N,P);
end

[THETA,fval,exitflag] = fmincon(fun,THmax,H,b,[],[],LB,UB,[],ops);

THETA*3600

if exitflag>=1
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
    Constraints_LWR_QFP_scr

    kk=2;

    while (norm(QH{kk}-QH{kk-1}))> 200 && exitflag>=1 && kk<6
        
        [THETA,fval,exitflag] = fmincon(fun,THmax,H,b,[],[],LB,UB,[],ops);

        THETA*3600;

        disp('end time is')
        sum(C)*3600
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
        
        Delay_params_LWR_QD_scr
        Constraints_LWR_QFP_scr

        disp('qa is: ')
        qa
        
    end
    
    norm(QH{kk}-QH{kk-1})
    [THETA,fval,exitflag] = fmincon(fun,THmax,H,b,[],[],LB,UB,[],ops);
    THETA*3600
   
end

%% TOTAL DELAY AND RESIDUAL QUEUES FOR EACH CYCLE

delta1 = zeros(P,N);
delta2 = zeros(P,N);
delta3 = zeros(P,N);
delta4 = zeros(P,N);

delay_pk = zeros(P,N); % delay of each phase in each cycle
Mixed_obj = zeros(P,N); % the mixed delay and spillback probability objective function
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
        
        Eps = [1,1];
        
        Mixed_obj(p,k) = omega1*delay_pk(p,k)+ 0.5*omega2*(Eps(p)/(Alphas(p)*Del(p)-delta2(p,k))+Eps(p)/(Alphas(p)*Del(p)-delta4(p,k)));

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
    info_struct_QF(N).delta1 = delta1;
    info_struct_QF(N).delta2 = delta2;
    info_struct_QF(N).delta3 = delta3;
    info_struct_QF(N).delta4 = delta4;
    info_struct_QF(N).C = C;
    info_struct_QF(N).TT = TT;
    info_struct_QF(N).qa = qa;
    info_struct_QF(N).delay = delay_pk;
    info_struct_QF(N).through = TH_pk;
    info_struct_QF(N).Mixed_obj = Mixed_obj;
    info_struct_QF(N).Xj = Xj;
    info_struct_QF(N).Theta = THETA;
    info_struct_QF(N).exitflag = exitflag;
    info_struct_QF(N).CritTime = CritTime;
    info_struct_QF(N).AN = AN;
    info_struct_QF(N).BN = BN;

else
    info_struct_QF(N).exitflag = exitflag;
end