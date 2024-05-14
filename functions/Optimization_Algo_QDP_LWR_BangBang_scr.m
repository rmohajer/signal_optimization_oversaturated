% This code performs the optimization algorithm for queue discharging
% period and using the LWR theory and bang bang strategy

%% ARRIVAL FLOWS MAPPING TO FUTURISTIC FLOWS


CH{1} = ones(N,1)*cmax; % history of cycle times

C = CH{1};

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

Constraints_LWR_QD_BB_scr

%% OPTIMIZATION OF THE TOTAL DELAY


LB = THmin;
UB = THmax;
% LB=[];
% UB=[];
ops = optimoptions('fmincon','Algorithm','sqp'); % active-set OR interior-point-convex

% [x_solN,fval,exitflag] = quadprog(2*AN,BN,H,b,Heq,beq,LB,UB,THmin,ops)

fun = @(x) x'*AN*x+BN'*x+CN;
% fun = @(x)BN'*x+CN;
[THETAbb,fval,exitflag] = fmincon(fun,THmin,H,b,Heq,beq,LB,UB,[],ops);
exitflag


CH{2} = zeros(N,1);

for j=1:N
    CH{2}(j) = sum(THETAbb((j-1)*P+1:j*P))+L;
end

C = CH{2};

Nf = N;
FlowMapping_scr
QH{2} = qa;

Delay_params_LWR_QD_scr
Constraints_LWR_QD_BB_scr

kk=2;

while (norm(CH{kk}-CH{kk-1}))>5/3600*N && exitflag==1 && kk<6

    Nf = N;
    FlowMapping_scr
    QH{kk} = qa;

    Delay_params_LWR_QD_scr
    Constraints_LWR_QD_BB_scr
    
    [THETAbb,fval,exitflag] = fmincon(fun,THmin,H,b,Heq,beq,LB,UB,[],ops);
    
    kk=kk+1;
    CH{kk} = zeros(N,1);

    for j=1:N
        CH{kk}(j) = sum(THETAbb((j-1)*P+1:j*P))+L;
    end

    C = CH{kk};
    
end


% [THETA,fval,exitflag] = quadprogIP(2*AN,BN,H,b,Heq,beq,LB,UB)
%  [THETA,fval,exitflag] = cplexqp([],BN,H,b,Heq,beq,LB,UB)
%  THETA*3600
%% TOTAL DELAY AND RESIDUAL QUEUES FOR EACH CYCLE

if exitflag==1
    delta1 = zeros(P,N);
    delta2 = zeros(P,N);
    delta3 = zeros(P,N);
    delta4 = zeros(P,N);

    delay_pk = zeros(P,N); % delay of each phase in each cycle
    Xj = zeros(P,N); % back of the queue in each cycle and each phase

    CritTime  = zeros(P,3*N);


    rownum=1;
    for k=1:N

        for p=1:P
            delta1(p,1) = deltBB_0(p);

            % calculation of r_tilde_p,1(k) and r_tilde_p,2(k)
            r_tilpk1 = 0;
            r_tilpk2 = 0;

            if p>1
                for j=1:(p-1)
                    r_tilpk1 = r_tilpk1+THETAbb((k-1)*P+j)+l(j);
                end % for j
            end

            if p<P
                for j=(p+1):P
                    r_tilpk2 = r_tilpk2+THETAbb((k-1)*P+j)+l(j);
                end
            end
            r_tilpk2 = r_tilpk2+l(p);



            if k>1
                delta1(p,k) = delta4(p,k-1);
            end

            delta2(p,k) = delta1(p,k)+Gamma(p,k)*r_tilpk1;
            delta3(p,k) = delta2(p,k)-qc(p)/kj(p)*THETAbb((k-1)*P+p);
            delta4(p,k) = delta3(p,k)+Gamma(p,k)*r_tilpk2;

            delay_pk(p,k) = kj(p)*(delta1(p,k)+0.5*(delta2(p,k)-delta1(p,k)))*r_tilpk1...
                +kj(p)*(delta3(p,k)+0.5*(delta4(p,k)-delta3(p,k)))*r_tilpk2;

            Xj(p,k) = max([delta2(p,k),delta4(p,k)]); 
            rownum = rownum+1;

            % calculation of the critical points
             if k==1
                CritTime(p,(k-1)*4+1) = Time;
                CritTime(p,(k-1)*4+2) = CritTime(p,(k-1)*4+1)+r_tilpk1;
                CritTime(p,(k-1)*4+3) = CritTime(p,(k-1)*4+2)+THETAbb((k-1)*P+p);
                CritTime(p,(k-1)*4+4) = CritTime(p,(k-1)*4+3)+r_tilpk2;
            else
                CritTime(p,(k-1)*4+1) = CritTime(p,(k-1)*4);
                CritTime(p,(k-1)*4+2) = CritTime(p,(k-1)*4+1)+r_tilpk1;
                CritTime(p,(k-1)*4+3) = CritTime(p,(k-1)*4+2)+THETAbb((k-1)*P+p);
                CritTime(p,(k-1)*4+4) = CritTime(p,(k-1)*4+3)+r_tilpk2;
            end


        end % for p
    end % for k
end

disp('delta3:')
delta3

% if bbloop==0
%     CritTime = [zeros(2,1),CritTime];
% end

TT = sum(C)*60; % total time spent

if exitflag==1
    BBexitflag = 1;
    % saving the output structure
    info_structBB.delta1 = [info_structBB.delta1,delta1];
    info_structBB.delta2 = [info_structBB.delta2,delta2];
    info_structBB.delta3 = [info_structBB.delta3,delta3];
    info_structBB.delta4 = [info_structBB.delta4,delta4];
    info_structBB.C = [info_structBB.C;C];
    info_structBB.TT = [info_structBB.TT,TT];
    info_structBB.qa = [info_structBB.qa,qa];
    info_structBB.delay = [info_structBB.delay,delay_pk];
    info_structBB.Xj = [info_structBB.Xj,Xj];
    info_structBB.Theta = [info_structBB.Theta;THETAbb];
    info_structBB.Eta = Eta;
    info_structBB.CritTime = [info_structBB.CritTime,CritTime];
    info_structBB.Nopt = [info_structBB.Nopt,N];
else
    BBexitflag=0 ;
end