% This code performs the optimization algorithm for queue discharging
% period and using the LWR theory and discharging the queues in a single
% cycle

%% ARRIVAL FLOWS MAPPING TO FUTURISTIC FLOWS
cmax = 200/3600; % h
CH{1} = ones(1,1)*cmax; % history of cycle times

C = CH{1};
Time = 0;
Nf = 1;
FlowMapping_scr
QH{1} = qa;

delt1_0=60/kj1; 
delt2_0=40/kj2;

delt_0 = [delt1_0;delt2_0];


%% CALCULATION OF THE BASIC MATRICES

% the following script calculates the parameters of the delay function for
% optimization.

A1 = zeros(P);
Ap = {};
for p=1:P
    Ap{p} = zeros(P);
    for i=1:P
        for j=1:P

            
            if i~=p && j~=p
                Ap{p}(i,j)=0.5*Gamma(p,1)*kj(p);
            elseif i==j && i==p
                Ap{p}(i,j)=0;
            else % if i=p or j=p
                if i>p || j>p
                    Ap{p}(i,j) = 0;
                end
            end
            
        end % for j
    end % for i
end % for p
A1 = plus(Ap{:});

B1 = zeros(P,1);

Bp = {};
for p = 1:P
    for i = 1:P
        L1=0;
        if p>1
            for j=1:p-1
                L1 = L1+l(j);
            end
        end

        L2 = 0;
        for j=p:P
            L2=L2+l(j);
        end
        
        if i<p
            Bp{p}(i,1) = (delt_0(p)+Gamma(p,1)*L1)*kj(p);
        elseif i>p
            Bp{p}(i,1) = (Gamma(p,1)*L2)*kj(p);
        else 
            Bp{p}(i,1) = 0;
        end
    end %for i
end % for p
B1 = plus(Bp{:});


%% CONSTRAINTS 

% the following script calculates the parameters of the constraints for
% optimization.

H4 = zeros(P,P); % for equality constraint
b4 = zeros(P,1);

% costruction of H4
for p = 1:P
    Hp4 = zeros(1,P);
    
    Htil_p4 = zeros(1,P);

    for jp=1:P
        if jp<p
            Htil_p4(jp) = Gamma(p,1);
        elseif jp==p
            Htil_p4(jp) = -qc(p)/kj(p);
        else
            Htil_p4(jp) = 0;
        end % if j'
    end % for j'

    Hp4(1,1:P) = Htil_p4;
    
    H4(p,:) = Hp4;
end % for p


% Construction of b4
for p=1:P
    
    bpk4 = -delt_0(p);
    
    if p>1
        for j=1:p-1
            bpk4 = bpk4- Gamma(p,1)*l(j);
        end % for j
    end % if p
    b4(p,1) = bpk4;
end % for p

THmin = gmin;
THmax = gmax;

% construction of H5 and b5
H5 = -ones(1,P);% for maximum cycle length
b5 = L-cmax;


% construction of H6 and b6
H6 = zeros(P); % for undersaturation constraint
b6 = zeros(P,1);

rownum = 1;

for p=1:P
    Hp6 = zeros(1,P);
    
    for jp=1:P
        if jp~=p
            Hp6(jp)=1-Eta(p,1)-1/(P-1);
        else
            Hp6(jp)=1;
        end
    end
    
    H6(rownum,:) = Hp6;
    rownum = rownum+1;
end % for k


rownum = 1;
for p=1:P
    b6(rownum,1) = -(1-Eta(p,1)-1/(P-1))*L-(2-P)/(P-1)*l(p);
    rownum = rownum+1;
end

% construction of H7 and b7
H7 = zeros((P-1),P); % for spillback avoidance constraint
b7 = zeros((P-1),1);

rownum = 1;

for p=2:P
    Hp7 = zeros(1,P);

    for jp=1:P
        if jp<p
            Hp7(jp) = -Gamma(p,1);
        else
            Hp7(jp) = 0;
        end % if jp
    end % for jp

    H7(rownum,:) = Hp7;
    rownum = rownum+1;
end % for k


rownum = 1;
 
for p=2:P
    b7(rownum,1) = -Del(p)+delt_0(p);

    for j=1:(p-1)
         b7(rownum,1) = b7(rownum,1)+Gamma(p,1)*l(j);

    end % for j
    rownum = rownum+1;
end
   

% construction of H8 and b8

H8 = zeros((P-1),P); % for spillback avoidance constraint
b8 = zeros((P-1),1);

rownum = 1;

for p=1:(P-1)
    Hp8 = zeros(1,P);

    for jp=1:P
        if jp~=p
            Hp8(jp) = -Gamma(p,1);
        else
            Hp8(jp) = qc(p)/kj(p);
        end % if jp
    end % for jp

    H8(rownum,:) = Hp8;
    rownum = rownum+1;
end % for k

rownum = 1;
    
for p=1:(P-1)
    b8(rownum,1) = -Del(p)+delt_0(p);
    b8(rownum,1) = b8(rownum,1)+Gamma(p,1)*L;

    rownum = rownum+1;
end

% H = -[-H4;H5;H6;H7;H8];
% b = -[-b4;b5;b6;b7;b8];

H = -[-H4;H6;H7;H8];
b = -[-b4;b6;b7;b8];

%% OPTIMIZATION OF THE TOTAL DELAY using the generic method


LB = THmin;
UB = THmax;
ops = optimoptions('fmincon','Algorithm','sqp'); % active-set OR interior-point-convex

% [x_solN,fval,exitflag] = quadprog(2*AN,BN,H,b,Heq,beq,LB,UB,THmin,ops)

fun = @(x) x'*A1*x+B1'*x;
% fun = @(x) B1'*x;
[THETA,fval,exitflag] = fmincon(fun,THmin,H,b,[],[],LB,[],[],ops)

C = sum(THETA)+L;
%% OPTIMIZATION OF THE TOTAL DELAY using the analytical method for a two-phase approach

Gtilde = [(Eta(1,1)/(1-Eta(1,1)*Eta(2,1))-1)*L; (Eta(1,1)*Eta(2,1)/(1-Eta(1,1)*Eta(2,1)))*L];

Gstar = [kj(1)/qc(1)*delt_0(1); Eta(2,1)*(kj(1)/qc(1)*delt_0(1)+l(1))+kj(1)/qc(2,1)*delt_0(2)];

Gopt = zeros(2,1);
if -H6*Gstar<=-b6
    Gopt = Gstar
    disp('step i. is passed')
elseif Gstar(1)>=Gtilde(1) && kj(2)/qc(2)*delt_0(2)+Eta(2,1)*l(1)>Eta(2,1)*L
    Gopt(2) = 1/(1-Eta(1,1)*Eta(2,1))*(Eta(1,1)*Eta(2,1)*L+Eta(2,1)*l(1)+kj(2)/qc(2)*delt_0(2));
    Gopt(1) = Eta(1,1)*(Gopt(2)+L);
    Gopt
    disp('step ii. is passed')
elseif Gstar(1)<=Gtilde(1) && kj(2)/qc(2)*delt_0(2)+Eta(2,1)*l(1)<= Eta(2,1)*L
    Gopt = Gtilde
    disp('step iii. is passed')
elseif Gstar(1)>Gtilde(1) && delt_0(2)+Eta(2,1)*l(1)<= Eta(2,1)*L
    Gopt = [kj(1)/qc(1)*delt_0(1); Eta(2,1)*(kj(1)/qc(1)*delt_0(1)+L)]
    disp('step iv. is passed')
else
    disp('It looks that there is no solution')
end


-H5*Gopt<=-b5

% THETA = Gopt;

%% TOTAL DELAY AND RESIDUAL QUEUES FOR EACH CYCLE

delta1 = zeros(P,1);
delta2 = zeros(P,1);
delta3 = zeros(P,1);
delta4 = zeros(P,1);

delay_pk = zeros(P,1); % delay of each phase in each cycle
Xj = zeros(P,1); % back of the queue in each cycle and each phase

CritTime  = zeros(P,3);


rownum=1;

for p=1:P
    delta1(p,1) = delt_0(p);

    % calculation of r_tilde_p,1(k) and r_tilde_p,2(k)
    r_tilpk1 = 0;
    r_tilpk2 = 0;

    if p>1
        for j=1:(p-1)
            r_tilpk1 = r_tilpk1+THETA(j)+l(j);
        end % for j
    end

    if p<P
        for j=(p+1):P
            r_tilpk2 = r_tilpk2+THETA(j)+l(j);
        end
    end
    r_tilpk2 = r_tilpk2+l(p);


    delta2(p,1) = delta1(p,1)+Gamma(p,1)*r_tilpk1;
    delta3(p,1) = delta2(p,1)-qc(p)/kj(p)*THETA(p);
    delta4(p,1) = Gamma(p,1)*r_tilpk2;

    delay_pk(p,1) = kj(p)*(delta1(p,1)+0.5*(delta2(p,1)-delta1(p,1)))*r_tilpk1...
        +kj(p)*(0.5*delta4(p,1))*r_tilpk2;

    Xj(p,k) = max([delta2(p,1),delta4(p,1)]); 
    rownum = rownum+1;

    % calculation of the critical points
    CritTime(p,1) = r_tilpk1;
    CritTime(p,2) = CritTime(p,1)+THETA(p);
    CritTime(p,3) = CritTime(p,2)+r_tilpk2;


end % for p


TT = sum(C)*60; % total time spent in [min]

if exitflag==1
    % saving the output structure
    info_struct(1).delta1 = delta1;
    info_struct(1).delta2 = delta2;
    info_struct(1).delta3 = delta3;
    info_struct(1).delta4 = delta4;
    info_struct(1).C = C;
    info_struct(1).TT = TT;
    info_struct(1).qa = qa;
    info_struct(1).delay = delay_pk;
    info_struct(1).Xj = Xj;
    info_struct(1).Theta = THETA;
    info_struct(1).exitflag = exitflag;
    info_struct(1).CritTime = CritTime;
    info_struct(1).AN = A1;
    info_struct(1).BN = B1;
    info_struct(1).Eta = Eta;

else
    info_struct(1).exitflag = exitflag;
end