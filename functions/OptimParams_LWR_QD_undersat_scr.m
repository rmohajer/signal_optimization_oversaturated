% this code calculates the parameters of the optimization of delay during undersaturation period
% The delay parameters and constraint parameters are derived.







%% Delay parameters

ANt_u = zeros(Nt*P);

for n1 = 1:Nt
    for n2 = 1:Nt
        Ap = {};
        for p=1:P
            Ap{p} = zeros(P);
            for i=1:P
                for j=1:P
                    if n1==n2
                        if i~=p && j~=p
                            Ap{p}(i,j)=0.5*Gamma(p,n1)*kj(p);
                        else
                            Ap{p}(i,j)=0;
                        end

                    else % if n1~=n2
                        if n1 == (n2-1)

                            if i>p && j<p
                                Ap{p}(i,j) = 0.5*Gamma(p,n1)*kj(p);
                            else % i<=p or j>=p
                                Ap{p}(i,j)=0;
                            end
                        end % if n1==n2-1
                    end
                end % for j
            end % for i
        end % for p
        Atp = plus(Ap{:});
        ANt_u((n1-1)*P+1:n1*P,(n2-1)*P+1:n2*P) = Atp;
    end % for n2
    
end % for n1

%%
BNt_u = zeros(Nt*P,1);

for k = 1:Nt
    Bp = {};
    for p = 1:P
        for j = 1:P
            Bp{p}(j) = 0;
            if j<p
                if p>1
                    for i=1:(p-1)
                        Bp{p}(j) = Bp{p}(j)+(Gamma(p,k)*l(i))*kj(p);
                    end
                end
                if k>1
                    for i=p:P
                        Bp{p}(j) = Bp{p}(j)+(Gamma(p,k-1)*l(i))*kj(p);
                    end
                end
            elseif j>p
                if k<Nt
                    Bp{p}(j)=Gamma(p,k)*L*kj(p);
                else % if k==Nt
                    for i=p:P
                        Bp{p}(j) = Bp{p}(j)+Gamma(p,k)*l(i)*kj(p);
                    end
                end
                
            else 
                Bp{p}(j) = 0;
            end
        end %for i
    end % for p
    Btp = plus(Bp{:});
    BNt_u((k-1)*P+1:k*P) = Btp;
end % for k

%% Constraint parameters

THmin_u = gmin;
for n=1:(Nt-1)
    THmin_u = [THmin_u;gmin];
end

THmax_u = gmax;
for n=1:(Nt-1)
    THmax_u = [THmax_u;gmax];
end


% construction of H5 and b5

H5_u = zeros(Nt, Nt*P); % for maximum cycle length
b5_u = zeros(Nt,1);

for k=1:Nt
    H5_u(k,:) = [zeros(1,(k-1)*P),-ones(1,P),zeros(1,(Nt-k)*P)];
    b5_u(k,1) = L-cmax;
end

% construction of H6 and b6

H6_u = zeros(Nt*P); % for undersaturation constraint
b6_u = zeros(Nt*P,1);

rownum = 1;
for k=1:Nt
    for p=1:P
        Hp6 = zeros(1,Nt*P);
        for kp=1:Nt
            H6pk = zeros(1,P);
            if kp==k
                for jp=1:P
                    if jp~=p
                        H6pk(jp)=1-Eta(p,k)-1/(P-1);
                    else
                        H6pk(jp)=1;
                    end
                end
                Hp6(1,(kp-1)*P+1:kp*P) = H6pk;
            else
                Hp6(1,(kp-1)*P+1:kp*P) = zeros(1,P);
            end
            
        end
        H6_u(rownum,:) = Hp6;
        rownum = rownum+1;
    end % for k
end % for p

rownum = 1;
for k=1:Nt
    for p=1:P
        b6_u(rownum,1) = -(1-Eta(p,k)-1/(P-1))*L-(2-P)/(P-1)*l(p);
        rownum = rownum+1;
    end
end

H_u = -[H5_u;H6_u];
b_u = -[b5_u;b6_u];

% H_u = -[H6_u];
% b_u = -[b6_u];

%% equality constraint \delta_{p,3}(k)=0

H1_u = zeros((Nt-1)*P,Nt*P);
b1_u = zeros((Nt-1)*P,1);

% costruction of H1
rownum = 1;
for k = 1:Nt
    for p = 1:P
        Hpk1 = zeros(1,Nt*P);
        for kp = 1:Nt
            Htil_pk1 = zeros(1,P);
            if kp<k
                for jp = 1:P
                    if jp~=p
                        Htil_pk1(jp) = Gamma(p,kp);
                    else
                        Htil_pk1(jp) = -qc(p)/kj(p);
                    end % if jp
                end % for jp
                
            elseif kp==k
                for jp=1:P
                    if jp<p
                        Htil_pk1(jp) = Gamma(p,kp);
                    elseif jp==p
                        Htil_pk1(jp) = -qc(p)/kj(p);
                    else
                        Htil_pk1(jp) = 0;
                    end % if j'
                end % for j'
                
            else % if k'>k
                Htil_pk1 = zeros(1,P);
            end % if k'
            Hpk1(1,(kp-1)*P+1:kp*P) = Htil_pk1;
        end % for k'
        H1_u(rownum,:) = Hpk1;
        rownum = rownum+1;
    end % for k
end % for p
                    

% Construction of b1

rownum = 1;
for k=1:Nt
    for p=1:P
        bpk1 = -info_struct(Np).delta4(p,end);
        if k>1
            for m=1:k-1
                bpk1 = bpk1-Gamma(p,m)*L;
            end % for m
        end %if k
        
        if p>1
            for j=1:p-1
                bpk1 = bpk1-Gamma(p,k)*l(j);
            end % for j
        end % if p
        b1_u(rownum,1) = bpk1;
        rownum = rownum+1;
    end % for k
end % for p

Heq_u = H1_u;
beq_u = b1_u;