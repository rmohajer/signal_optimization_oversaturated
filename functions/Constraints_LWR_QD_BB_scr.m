
H1_bb = zeros((N)*(P),N*P); % res queues must be larger than zero
H1t_bb = zeros((N)*(P-1),N*P); %  res queue length stays constant equality constraint 
% H3_bb = zeros((N-1)*(P-1),N*P); % for non-decreasing res queue length equality constraint
H3p_bb = zeros((N-1),N*P); % for non-decreasing res queue length at phase p_bb
H4_bb = zeros(1,N*P); % for equality constraint at phase p_bb

b1_bb = zeros(N,1);
b1t_bb = zeros(N,1);
% b3_bb = zeros((N-1)*(P-1),1);
b3p_bb = zeros((N-1),1);
b4_bb = zeros(1,1);


% costruction of H1_bb

rownum = 1;
for k = 1:N
    for p = P_larger
        Hpk1 = zeros(1,N*P);
        for kp = 1:N
            Htil_pk1 = zeros(1,P);
            if kp<k
                for jp = 1:P
                    if jp~=p
                        Htil_pk1(jp) = Gamma(p,k0+kp);
                    else
                        Htil_pk1(jp) = -qc(p)/kj(p);
                    end % if jp
                end % for jp
                
            elseif kp==k
                for jp=1:P
                    if jp<p
                        Htil_pk1(jp) = Gamma(p,k0+kp);
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
        H1_bb(rownum,:) = Hpk1;
        rownum = rownum+1;
    end % for k
end % for p
                    

% Construction of b1_bb

rownum = 1;
for k=1:N
    for p=1:P
        bpk1 = -delt_0(p);
        if k>1
            for m=1:k-1
                bpk1 = bpk1-Gamma(p,k0+m)*L;
            end % for m
        end %if k
        
        if p>1
            for j=1:p-1
                bpk1 = bpk1-Gamma(p,k0+k)*l(j);
            end % for j
        end % if p
        b1_bb(rownum,1) = bpk1;
        rownum = rownum+1;
    end % for k
end % for p




% costruction of H1t_bb

p_select = 1:P;
p_select(p_select==p_bb) = [];

rownum = 1;
for k = 1:N
    for p=p_select
        Hpk1 = zeros(1,N*P);
        for kp = 1:N
            Htil_pk1 = zeros(1,P);
            if kp<k
                for jp = 1:P
                    if jp~=p
                        Htil_pk1(jp) = Gamma(p,k0+kp);
                    else
                        Htil_pk1(jp) = -qc(p)/kj(p);
                    end % if jp
                end % for jp

            elseif kp==k
                for jp=1:P
                    if jp<p
                        Htil_pk1(jp) = Gamma(p,k0+kp);
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
        H1t_bb(rownum,:) = Hpk1;
        rownum = rownum+1;
    end % for p
end % for k
                    

% Construction of b1t_bb

rownum = 1;
for k=1:N
    for p=p_select
        
        bpk1 = -deltBB_0(p)+Res_que0(p);
        if k>1
            for m=1:k-1
                bpk1 = bpk1-Gamma(p,k0+m)*L;
            end % for m
        end %if k

        if p>1
            for j=1:p-1
                bpk1 = bpk1-Gamma(p,k0+k)*l(j);
            end % for j
        end % if p
        b1t_bb(rownum,1) = bpk1;
        rownum = rownum+1;
    end % for p
    
end % for k





% costruction of H3p_bb

rownum = 1;
for k = 1:N-1
    
    Hpk3 = zeros(1,N*P);
    for kp = 1:N
        Htil_pk3 = zeros(1,P);
        if kp==k
            for jp = 1:P
                if jp<=p_bb
                    Htil_pk3(jp) = 0;
                else
                    Htil_pk3(jp) = -Gamma(p_bb,k0+kp);
                end % if jp
            end % for jp

        elseif kp==k+1
            for jp=1:P
                if jp<p_bb
                    Htil_pk3(jp) = -Gamma(p_bb,k0+kp);
                elseif jp==p_bb
                    Htil_pk3(jp) = qc(p_bb)/kj(p_bb);
                else
                    Htil_pk3(jp) = 0;
                end % if j'
            end % for j'

        else % if k'>k
            Htil_pk3 = zeros(1,P);
        end % if k'
        Hpk3(1,(kp-1)*P+1:kp*P) = Htil_pk3;
    end % for k'
    H3p_bb(rownum,:) = Hpk3;
    rownum = rownum+1;
   
end % for k


% Construction of b3p_bb

rownum = 1;
for k=1:N-1
    
    bpk3 = 0;
    for j=p_bb:P
        bpk3 = bpk3+Gamma(p_bb,k0+k)*l(j);
    end % for j

    if p>1
        for j=1:p-1
            bpk3 = bpk3+ Gamma(p_bb,k0+k+1)*l(j);
        end % for j
    end % if p
    b3p_bb(rownum,1) = bpk3;
    rownum = rownum+1;
    
end % for k


% costruction of H4_bb


Hp4 = zeros(1,N*P);
for kp = 1:N
    Htil_p4 = zeros(1,P);
    if kp<N
        for jp = 1:P
            if jp~=p_bb
                Htil_p4(jp) = Gamma(p_bb,k0+kp);
            else
                Htil_p4(jp) = -qc(p_bb)/kj(p_bb);
            end % if jp
        end % for jp

    else
        for jp=1:P
            if jp<p_bb
                Htil_p4(jp) = Gamma(p_bb,k0+kp);
            elseif jp==p_bb
                Htil_p4(jp) = -qc(p_bb)/kj(p_bb);
            else
                Htil_p4(jp) = 0;
            end % if j'
        end % for j'
    end % if k'
    Hp4(1,(kp-1)*P+1:kp*P) = Htil_p4;
end % for k'
H4_bb = Hp4;



% Construction of b4


    
bpk4 = -deltBB_0(p_bb);
for m=1:(N-1)
    bpk4 = bpk4-Gamma(p_bb,k0+m)*L;
end % for m

if p_bb>1
    for j=1:p_bb-1
        bpk4 = bpk4- Gamma(p_bb,k0+N)*l(j);
    end % for j
end % if p
b4_bb = bpk4;


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
                        H6pk(jp)=1-Eta(p,k0+k)-1/(P-1);
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
        b6(rownum,1) = -(1-Eta(p,k0+k)-1/(P-1))*L-(2-P)/(P-1)*l(p);
        rownum = rownum+1;
    end
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
                        H7pk(jp) = -Gamma(p,k0+kp);
                    else
                        H7pk(jp) = qc(p)/kj(p);
                    end % if jp
                end % for jp
            elseif kp==k
                for jp=1:P
                    if jp<p
                        H7pk(jp) = -Gamma(p,k0+kp);
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
        b7(rownum,1) = -Del(p)+deltBB_0(p);
        for m=1:(k-1)
            b7(rownum,1) = b7(rownum,1)+Gamma(p,k0+m)*L;
        end % for m 

        for j=1:(p-1)
             b7(rownum,1) = b7(rownum,1)+Gamma(p,k0+k)*l(j);
             
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
                        H8pk(jp) = -Gamma(p,k0+kp);
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
        b8(rownum,1) = -Del(p)+deltBB_0(p);
        for m=1:k
            b8(rownum,1) = b8(rownum,1)+Gamma(p,k0+m)*L;
        end % for m 
        rownum = rownum+1;
    end
   
end



% H = -[H3p_bb;H5;H7;H8;H6];
% b = -[b3p_bb;b5;b7;b8;b6];

H = -[H3p_bb;H5;H7;H8];
b = -[b3p_bb;b5;b7;b8];

Heq = [H4_bb;H1t_bb];
beq = [b4_bb;b1t_bb];

% Heq = [H1t_bb];
% beq = [b1t_bb];