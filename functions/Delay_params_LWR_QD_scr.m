
% This script calculates the Matrices AN, BN, and CN that result in the
% delay function from the LWR theory based model

AN = zeros(N*P);

for n1 = 1:N
    for n2 = 1:N
        Ap = {};
        
        for p=1:P
            Ap{p} = zeros(P);
            if n1~=N && n2~=N
                for i=1:P
                    for j=1:P

                        if n1==n2
                            if i~=p && j~=p
                                Ap{p}(i,j)=0.5*Gamma(p,n1)*kj(p);
                            elseif i==j && i==p
                                Ap{p}(i,j)=0;
                            else % if i=p or j=p
                                if i>p || j>p
                                    Ap{p}(i,j) = -0.5*qc(p)/kj(p)*kj(p);
                                end
                            end

                        else % if n1~=n2
                            n = min(n1,n2);
                            if n1<n2
                                
                                if i~=p && j~=p
                                    Ap{p}(i,j) = 0.5*Gamma(p,n)*kj(p);
                                elseif i==p && j~=p
                                    Ap{p}(i,j) = -0.5*qc(p)/kj(p)*kj(p);
                                else %  j=p
                                    Ap{p}(i,j)=0;
                                end
                            else % if n1>n2
                                if i~=p && j~=p
                                    Ap{p}(i,j) = 0.5*Gamma(p,n)*kj(p);
                                elseif j==p && i~=p
                                    Ap{p}(i,j) = -0.5*qc(p)/kj(p)*kj(p);
                                else %  i=p
                                    Ap{p}(i,j)=0;
                                end
                            end
                        end
                    end % for j
                end % for i
            else % if ni=N
                for i=1:P
                    for j=1:P
                        
                        if n1==n2 % =N
                            if i~=p && j~=p
                                Ap{p}(i,j)=0.5*Gamma(p,n1)*kj(p);
                            else % if i==p or j==p
                                Ap{p}(i,j)=0;
                            end

                        else % if n1~=n2
                            n = min(n1,n2);
                            if n1<n2
                                if i~=p && j<p
                                    Ap{p}(i,j)=0.5*Gamma(p,n)*kj(p);
                                elseif i==p && j<p
                                    Ap{p}(i,j)=-0.5*qc(p)/kj(p)*kj(p);
                                else % j>p
                                    Ap{p}(i,j)=0;
                                end
                            else % n1>n2
                               if j~=p && i<p
                                    Ap{p}(i,j)=0.5*Gamma(p,n)*kj(p);
                                elseif j==p && i<p
                                    Ap{p}(i,j)=-0.5*qc(p)/kj(p)*kj(p);
                                else % i>p
                                    Ap{p}(i,j)=0;
                                end 
                            end
                            
                        end
                    end % for j
                end % for i
                
            end % if ni~=N
        end % for p
        
        Atp = plus(Ap{:});
        AN((n1-1)*P+1:n1*P,(n2-1)*P+1:n2*P) = Atp;
    end % for n2
end % for n1

BN = zeros(N*P,1);

for k = 1:N
    Bp = {};
    if k<N
        for p = 1:P
            for i = 1:P

                if i~=p
                    L1 = 0;
                    if p>1
                        for j=1:p-1
                            L1=L1+l(j);
                        end % for j
                    end % if p>1
                    
                    if k>1
                        Bp{p}(i) = (delt_0(p)+sum(Gamma(p,1:k-1))*L+(N-k)*Gamma(p,k)*L+Gamma(p,k)*L1)*kj(p);
                    else
                        Bp{p}(i) = (delt_0(p)+(N-k)*Gamma(p,k)*L+Gamma(p,k)*L1)*kj(p);
                    end % if k>1
                else 
                    Bp{p}(i) = -(N-k)*qc(p)/kj(p)*L*kj(p);
                end
            end %for i
        end % for p
    else % if k=N
        for p = 1:P
            for i = 1:P
                L1 = 0;
                if p>1
                    for j=1:p-1
                        L1=L1+l(j);
                    end % for j
                end % if p>1
                
                L2 = 0;
                for j=p:P
                    L2=L2+l(j);
                end % for j

                if i<p
                    Bp{p}(i) = (delt_0(p)+sum(Gamma(p,1:k-1))*L+Gamma(p,k)*L1)*kj(p);
                elseif i>p
                    Bp{p}(i) = Gamma(p,k)*L2*kj(p);
                else
                    Bp{p}(i) =0;
                end
            end %for i
        end % for p
    end
    Btp = plus(Bp{:});
    BN((k-1)*P+1:k*P) = Btp;
end % for k

CN=0;
for k=1:N
    Cp = {};
    if k<N
        for p=1:P
            if k>1
                Cp{p}(k) = (delt_0(p)*L+sum(Gamma(p,1:(k-1)))*L^2+0.5*Gamma(p,k)*L^2)*kj(p);
            else
                Cp{p}(k) = (delt_0(p)*L+0.5*Gamma(p,k)*L^2)*kj(p);
            end
        end
    else % if k=N
        for p=1:P
            
            L1 = 0;
            if p>1
                for j=1:p-1
                    L1=L1+l(j);
                end % for j
            end % if p>1

            L2 = 0;
            for j=p:P
                L2=L2+l(j);
            end % for j
        
            Cp{p}(k) = ((delt_0(p)+sum(Gamma(p,1:(k-1)))*L)*L1+0.5*Gamma(p,k)*(L1^2+L2^2))*kj(p);
        end % for p
    end % if k<N
    Ctp = plus(Cp{:});
    
end
CN = sum(Ctp);  