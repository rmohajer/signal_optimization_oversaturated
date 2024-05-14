% this function calculates the spillback avoidance part of the objective
% function

function y = SP_fun(x,Del,delt_0,Gamma,qc,kj,l,N,P)

Eps = [1,1];
global Alphas


L = sum(l);

y = 0;

for k=1:N
    for p=1:P
        
        delta4_p = delt_0(p);
        
        for m=1:k
           delta4_p = delta4_p+ Gamma(p,m)*L-(Gamma(p,m)+qc(p)/kj(p))*x((m-1)*P+p);
           for j=1:P
               delta4_p = delta4_p+Gamma(p,m)*x((m-1)*P+j);
           end % for j
        end % for m
        
        delta2_p = delt_0(p);
        
        if k>1
            
            for m=1:k-1
                delta2_p = delta2_p+ Gamma(p,m)*L-(Gamma(p,m)+qc(p)/kj(p))*x((m-1)*P+p);
               for j=1:P
                   delta2_p = delta2_p+Gamma(p,m)*x((m-1)*P+j);
               end % for j
            end % for m
                
        end
        
        if p>1
            for j=1:p-1
                delta2_p = delta2_p+ Gamma(p,k)*(x((k-1)*P+j)+l(j));
            end
        end
        
        
        y = y+ Eps(p)/(Alphas(p)*Del(p)-delta2_p)+Eps(p)/(Alphas(p)*Del(p)-delta4_p);
        
    end % for p
end % for k


end

