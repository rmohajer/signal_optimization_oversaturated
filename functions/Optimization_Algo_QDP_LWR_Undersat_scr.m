
% This code performs the optimization algorithm for queue discharging
% period and using the LWR theory for undersaturated traffic mode
exitflag =1;
CH{1} = ones(Nt,1)*cmax; % history of cycle times

C = CH{1};

Nf = Nt;
FlowMapping_scr
QH{1} = qa;

OptimParams_LWR_QD_undersat_scr
LB = THmin_u;
UB = THmax_u;


fun_u = @(x) x'*ANt_u*x+BNt_u'*x;
% fun_u = @(x)repmat(qc,Nt,1)'*x;  % for throughput minimization
[THETA_u,fval,exitflag] = fmincon(fun_u,THmin_u,H_u,b_u,[],[],LB,UB,[],ops)
% [THETA_u,fval,exitflag,~] = quadprog(2*ANt_u,BNt_u,H_u,b_u,Heq_u,beq_u,LB,UB,THmin_u)

% THETA_u=Heq_u\beq_u;

CH{2} = zeros(Nt,1);

for j=1:Nt
    CH{2}(j) = sum(THETA_u((j-1)*P+1:j*P))+L;
end

C = CH{2};

Nf = Nt;
FlowMapping_scr
QH{2} = qa;

OptimParams_LWR_QD_undersat_scr

kk=2;

while (sum(CH{kk})-sum(CH{kk-1}))*60>1 

    [THETA_u,fval,exitflag] = fmincon(fun_u,THmin_u,H_u,b_u,[],[],LB,UB,[],ops)
%     [THETA_u,fval,exitflag,~] = quadprog(2*ANt_u,BNt_u,H_u,b_u,Heq_u,beq_u,LB,UB,THmin_u)

%     THETA_u=Heq_u\beq_u;
    kk=kk+1;
    CH{kk} = zeros(Nt,1);

    for j=1:Nt
        CH{kk}(j) = sum(THETA_u((j-1)*P+1:j*P))+L;
    end

    C = CH{kk};

    Nf = N;
    FlowMapping_scr
    QH{kk} = qa;

    OptimParams_LWR_QD_undersat_scr
    norm(QH{kk}-QH{kk-1})
end

[THETA_u,fval,exitflag] = fmincon(fun_u,THmin_u,H_u,b_u,[],[],LB,UB,[],ops)
% [THETA_u,fval,exitflag,~] = quadprog(2*ANt_u,BNt_u,H_u,b_u,Heq_u,beq_u,LB,UB,THmin_u)

% THETA_u=Heq_u\beq_u;

kk =kk+1;

CH{kk} = zeros(Nt,1);

for j=1:Nt
    CH{kk}(j) = sum(THETA_u((j-1)*P+1:j*P))+L;
end

C = CH{kk};
TT_u = sum(C)*60;


%%
if TT_u>= TTmax-info_struct(Np).TT
    flag_while=1;
    
    delta1 = zeros(P,Nt);
    delta3 = zeros(P,Nt);
    delta3v = zeros(P,Nt);
    delta2 = zeros(P,Nt);
    delta4 = zeros(P,Nt);

    delay_pk = zeros(P,Nt); % delay of each phase in each cycle
    TH_pk = zeros(P,Nt); % throughput of each phase in each cycle
    Xj = zeros(P,Nt); % back of the queue in each cycle and each phase

    CritTime  = zeros(P,3*Nt);


    rownum=1;
    for k=1:Nt

        for p=1:P
            delta1(p,1) = info_struct(Np).delta4(p,Np);

            % calculation of r_tilde_p,1(k) and r_tilde_p,2(k)
            r_tilpk1 = 0;
            r_tilpk2 = 0;
            

            if p>1
                for j=1:(p-1)
                    r_tilpk1 = r_tilpk1+THETA_u((k-1)*P+j)+l(j);
                end % for j
            end

            if p<P
                for j=(p+1):P
                    r_tilpk2 = r_tilpk2+THETA_u((k-1)*P+j)+l(j);
                end
            end
            r_tilpk2 = r_tilpk2+l(p);


            if k>1
                delta1(p,k) = delta4(p,k-1);
            end

            delta2(p,k) = delta1(p,k)+Gamma(p,k)*r_tilpk1;
            delta3v(p,k) = delta2(p,k)-qc(p)/kj(p)*THETA_u((k-1)*P+p);
            delta4(p,k) = Gamma(p,k)*r_tilpk2;

            delay_pk(p,k) = kj(p)*(delta1(p,k)+0.5*(delta2(p,k)-delta1(p,k)))*r_tilpk1...
                +kj(p)*(0.5*(delta4(p,k)))*r_tilpk2;
            
            TH_pk(p,k) = qc(p)*THETA_u((k-1)*P+p);


            Xj(p,k) = max([delta2(p,k),delta4(p,k)]); 
            rownum = rownum+1;

            % calculation of the critical points
            if k==1
                CritTime(p,(k-1)*4+1) = Time;
                CritTime(p,(k-1)*4+2) = CritTime(p,(k-1)*4+1)+r_tilpk1;
                CritTime(p,(k-1)*4+3) = CritTime(p,(k-1)*4+2)+THETA_u((k-1)*P+p);
                CritTime(p,(k-1)*4+4) = CritTime(p,(k-1)*4+3)+r_tilpk2;
            else
                CritTime(p,(k-1)*4+1) = CritTime(p,(k-1)*4);
                CritTime(p,(k-1)*4+2) = CritTime(p,(k-1)*4+1)+r_tilpk1;
                CritTime(p,(k-1)*4+3) = CritTime(p,(k-1)*4+2)+THETA_u((k-1)*P+p);
                CritTime(p,(k-1)*4+4) = CritTime(p,(k-1)*4+3)+r_tilpk2;
            end
            
            

        end % for p
    end % for k


    TT_u = sum(C)*60; % total time spent

    if exitflag==1
        % saving the output structure
        info_struct(Np).delta1_u = delta1;
        info_struct(Np).delta2_u = delta2;
        info_struct(Np).delta3_u = delta3;
        info_struct(Np).delta4_u = delta4;
        info_struct(Np).C_u = C;
        info_struct(Np).TT_u = TT_u;
        info_struct(Np).Nu = Nt;
        info_struct(Np).qa_u = qa;
        info_struct(Np).delay_u = delay_pk;
        info_struct(Np).through_u = TH_pk;
        info_struct(Np).Xj_u = Xj;
        info_struct(Np).Theta_u = THETA;
        info_struct(Np).exitflag_u = exitflag;
        info_struct(Np).CritTime_u = CritTime;
        info_struct(Np).AN_u = ANt_u;
        info_struct(Np).BN_u = BNt_u;

    else
        info_struct(Np).exitflag_u = exitflag;
    end


else
    Nt=Nt+1;
end