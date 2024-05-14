% This code calculates the constraint parameters and matrices using the LWR
% model for queue formation period (QFP)

H1 = zeros((N-1)*P,N*P); % non-negative res queue length

b1 = zeros((N-1)*P,1);

% costruction of H1
rownum = 1;
for k = 1:(N-1)
    for p = 1:P
        Hpk1 = zeros(1,N*P);
        for kp = 1:N
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
        H1(rownum,:) = Hpk1;
        rownum = rownum+1;
    end % for k
end % for p
                    

% Construction of b1

rownum = 1;
for k=1:N-1
    for p=1:P
        bpk1 = -delt_0(p);
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
        b1(rownum,1) = bpk1;
        rownum = rownum+1;
    end % for k
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


% construction of H7 and b7

H7 = zeros(N*(P-1),N*P); % for spillback avoidance constraint
b7 = zeros(N*(P-1),1);

rownum = 1;
for k=1:N
    for p=2:P
        Hp7 = zeros(1,N*P);
        
        for kp=1:N
            H7pk = zeros(1,P);
            if kp<k
                for jp=1:P
                    if jp~=p
                        H7pk(jp) = -Gamma(p,kp);
                    else
                        H7pk(jp) = qc(p)/kj(p);
                    end % if jp
                end % for jp
            elseif kp==k
                for jp=1:P
                    if jp<p
                        H7pk(jp) = -Gamma(p,kp);
                    else
                        H7pk(jp) = 0;
                    end % if jp
                end % for jp
            else % if kp>k
                H7pk = zeros(1,P);
            end % if kp
            Hp7(1,(kp-1)*P+1:kp*P) = H7pk;
        end % for kp
        
        
        H7(rownum,:) = Hp7;
        rownum = rownum+1;
    end % for k
end % for p

rownum = 1;
for k=1:N
    
    for p=2:P
        b7(rownum,1) = -Del(p)+delt_0(p);
        for m=1:(k-1)
            b7(rownum,1) = b7(rownum,1)+Gamma(p,m)*L;
        end % for m 

        for j=1:(p-1)
             b7(rownum,1) = b7(rownum,1)+Gamma(p,k)*l(j);
             
        end % for j
        rownum = rownum+1;
    end
   
end

% construction of H8 and b8

H8 = zeros(N*(P-1),N*P); % for spillback avoidance constraint
b8 = zeros(N*(P-1),1);

rownum = 1;
for k=1:N
    for p=1:(P-1)
        Hp8 = zeros(1,N*P);
        
        for kp=1:N
            H8pk = zeros(1,P);
            if kp<=k
                for jp=1:P
                    if jp~=p
                        H8pk(jp) = -Gamma(p,kp);
                    else
                        H8pk(jp) = qc(p)/kj(p);
                    end % if jp
                end % for jp
            else % if kp>k
                H8pk = zeros(1,P);
            end % if kp
            Hp8(1,(kp-1)*P+1:kp*P) = H8pk;
        end % for kp
        
        
        H8(rownum,:) = Hp8;
        rownum = rownum+1;
    end % for k
end % for p


rownum = 1;
for k=1:N
    
    for p=1:(P-1)
        b8(rownum,1) = -Del(p)+delt_0(p);
        for m=1:k
            b8(rownum,1) = b8(rownum,1)+Gamma(p,m)*L;
        end % for m 
        rownum = rownum+1;
    end
   
end


Heq = ones(1,N*P);
beq = T_total-(N-3)*L;

if flag_delay==2
    H = -[H1;H5;H7;H8;-Heq];
    b = -[b1;b5;b7;b8;-beq];
else

    H = -[H1;H5;H7;H8;Heq];
    b = -[b1;b5;b7;b8;beq];
end


% Heq = [];
% beq = [];
% 
% H9 = ones(1,N*P);
% b9 = T_total-N*L;
% 
% H = -[H1;H5;H7;H8;H9];
% b = -[b1;b5;b7;b8;b9];



