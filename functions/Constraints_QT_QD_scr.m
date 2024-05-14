% This code calculates the constraint parameters and matrices using the
% queuing theory based model for queue discharding (QD) period

H1q = zeros((N-1)*P,N*P); % non-negative res queue length
H2q = zeros((N-1)*P,N*P); % sum of res queue length larger than zero
H3q = zeros((N-1)*P,N*P); % for non-decreasing res queue length
H4q = zeros(P,N*P); % for equality constraint

b1q = zeros((N-1)*P,1);
b2q = zeros((N-1)*P,1);
b3q = zeros((N-1)*P,1);
b4q = zeros(P,1);

% costruction of H1
rownum = 1;
for k = 1:N-1
    for p = 1:P
        Hpk1 = zeros(1,N*P);
        for kp = 1:N
            Htil_pk1 = zeros(1,P);
            if kp<k
                for jp = 1:P
                    if jp~=p
                        Htil_pk1(jp) = qa(p,kp);
                    else
                        Htil_pk1(jp) = qa(p,kp)-qc(p);
                    end % if jp
                end % for jp
                
            elseif kp==k
                for jp=1:P
                    if jp<p
                        Htil_pk1(jp) = qa(p,kp);
                    elseif jp==p
                        Htil_pk1(jp) = qa(p,kp)-qc(p);
                    else
                        Htil_pk1(jp) = 0;
                    end % if j'
                end % for j'
                
            else % if k'>k
                Htil_pk1 = zeros(1,P);
            end % if k'
            Hpk1(1,(kp-1)*P+1:kp*P) = Htil_pk1;
        end % for k'
        H1q(rownum,:) = Hpk1;
        rownum = rownum+1;
    end % for k
end % for p
                    

% Construction of b1

rownum = 1;
for k=1:N-1
    for p=1:P
        bpk1 = -delt_0(p)*kj(p);
        if k>1
            for m=1:k-1
                bpk1 = bpk1-qa(p,m)*L;
            end % for m
        end %if k
        
        if p>1
            for j=1:p-1
                bpk1 = bpk1-qa(p,k)*l(j);
            end % for j
        end % if p
        b1q(rownum,1) = bpk1;
        rownum = rownum+1;
    end % for k
end % for p


% costruction of H2

rownum = 1;
for k = 1:N-1
    for p = 1:P
        Hpk2 = zeros(1,N*P);
        for kp = 1:N
            Htil_pk2 = zeros(1,P);
            if kp<k
                for jp = 1:P
                    if jp~=p
                        for r=1:P
                            Htil_pk2(jp) = Htil_pk2(jp)+qa(r,kp);
                        end % for r
                            
                    else
                        for r=1:P
                            Htil_pk2(jp) = Htil_pk2(jp)+qa(r,kp)-qc(r);
                        end % for r
                    end % if jp
                end % for jp
                
            elseif kp==k
                for jp=1:P
                    if jp<p
                        for r=1:P
                            Htil_pk2(jp) = Htil_pk2(jp)+qa(r,kp);
                        end % for r
                    elseif jp==p
                        for r=1:P
                            Htil_pk2(jp) = Htil_pk2(jp)+qa(r,kp)-qc(r);
                        end % for r
                    else
                        Htil_pk2(jp) = 0;
                    end % if j'
                end % for j'
                
            else % if k'>k
                Htil_pk2 = zeros(1,P);
            end % if k'
            Hpk2(1,(kp-1)*P+1:kp*P) = Htil_pk2;
        end % for k'
        H2q(rownum,:) = Hpk2;
        rownum = rownum+1;
    end % for k
end % for p


% Construction of b2

rownum = 1;
for k=1:N-1
    for p=1:P
        bpk2 = -P*delt_0(p)*kj(p);
        if k>1
            for r=1:P
                for m=1:k-1
                    bpk2 = bpk2-qa(r,m)*L;
                end % for m
            end % for r
        end %if k
        
        if p>1
            for r=1:P
                for j=1:p-1
                    bpk2 = bpk2-qa(r,k)*l(j);
                end % for j
            end % for r
        end % if p
        b2q(rownum,1) = bpk2;
        rownum = rownum+1;
    end % for k
end % for p

    
% costruction of H3

rownum = 1;
for k = 1:N-1
    for p = 1:P
        Hpk3 = zeros(1,N*P);
        for kp = 1:N
            Htil_pk3 = zeros(1,P);
            if kp==k
                for jp = 1:P
                    if jp<=p
                        Htil_pk3(jp) = 0;
                    else
                        Htil_pk3(jp) = -qa(p,kp);
                    end % if jp
                end % for jp
                
            elseif kp==k+1
                for jp=1:P
                    if jp<p
                        Htil_pk3(jp) = -qa(p,kp);
                    elseif jp==p
                        Htil_pk3(jp) = qc(p)-qa(p,kp);
                    else
                        Htil_pk3(jp) = 0;
                    end % if j'
                end % for j'
                
            else % if k'>k
                Htil_pk3 = zeros(1,P);
            end % if k'
            Hpk3(1,(kp-1)*P+1:kp*P) = Htil_pk3;
        end % for k'
        H3q(rownum,:) = Hpk3;
        rownum = rownum+1;
    end % for k
end % for p


% Construction of b3

rownum = 1;
for k=1:N-1
    for p=1:P
        bpk3 = 0;
        for j=p:P
            bpk3 = bpk3+qa(p,k)*l(j);
        end % for j
        
        if p>1
            for j=1:p-1
                bpk3 = bpk3+ qa(p,k+1)*l(j);
            end % for j
        end % if p
        b3q(rownum,1) = bpk3;
        rownum = rownum+1;
    end % for k
end % for p


% costruction of H4

for p = 1:P
    Hp4 = zeros(1,N*P);
    for kp = 1:N
        Htil_p4 = zeros(1,P);
        if kp<N
            for jp = 1:P
                if jp~=p
                    Htil_p4(jp) = qa(p,kp);
                else
                    Htil_p4(jp) = qa(p,kp)-qc(p);
                end % if jp
            end % for jp

        else
            for jp=1:P
                if jp<p
                    Htil_p4(jp) = qa(p,kp);
                elseif jp==p
                    Htil_p4(jp) = qa(p,kp)-qc(p);
                else
                    Htil_p4(jp) = 0;
                end % if j'
            end % for j'
        end % if k'
        Hp4(1,(kp-1)*P+1:kp*P) = Htil_p4;
    end % for k'
    H4q(p,:) = Hp4;
end % for p


% Construction of b4

for p=1:P
    
    bpk4 = -delt_0(p)*kj(p);
    for m=1:(N-1)
        bpk4 = bpk4-qa(p,m)*L;
    end % for m

    if p>1
        for j=1:p-1
            bpk4 = bpk4- qa(p,N)*l(j);
        end % for j
    end % if p
    b4q(p,1) = bpk4;
end % for p

THmin = gmin;
for n=1:(N-1)
    THmin = [THmin;gmin];
end

THmax = gmax;
for n=1:(N-1)
    THmax = [THmax;gmax];
end

% construction of H5 and b5

H5 = zeros(N, N*P); % for maximum cycle length
b5 = zeros(N,1);

for k=1:N
    H5(k,:) = [zeros(1,(k-1)*P),-ones(1,P),zeros(1,(N-k)*P)];
    b5(k,1) = L-cmax;
end

% construction of H6 and b6

H6 = zeros(N*P); % for undersaturation constraint
b6 = zeros(N*P,1);

rownum = 1;
for k=1:N
    for p=1:P
        Hp6 = zeros(1,N*P);
        for kp=1:N
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
        H6(rownum,:) = Hp6;
        rownum = rownum+1;
    end % for k
end % for p

rownum = 1;
for k=1:N
    for p=1:P
        b6(rownum,1) = -(1-Eta(p,k)-1/(P-1))*L-(2-P)/(P-1)*l(p);
        rownum = rownum+1;
    end
end







Hq = -[H1q;H2q;H3q;H5;H6];
bq = -[b1q;b2q;b3q;b5;b6];

% Hq = -[H1q;H2q];
% bq = -[b1q;b2q];

Heq_q = H4q;
beq_q = b4q;