% we need to map the flow to future flows 



qa = zeros(P,Nf);

for ppp=1:P
    Time_f = Time;
    n=floor(Time/ts)+1; % multiplier of ts showing that we need to switch to the next flows
    for k=1:Nf
         
        if Time_f > ts*n+ts/4
            n=n+1;
            qa(ppp,k) = Qa(ppp,n);   % max(Qa(ppp,n),Qa(ppp,n-1));
        else
            qa(ppp,k) = Qa(ppp,n);
        end
        Time_f = Time_f + C(k);
        disptime = Time_f*3600;
    end
end



% figure
% stairs(1:N,qa(1,:))
% hold on
% stairs(1:N,qa(2,:),'r')
% xlabel('cycle number')
% ylabel('q^a [veh/h]')
% legend('Approach 1','Approach 2')

% qa = [740*ones(1,N);520*ones(1,N)];

ka = qa./qc.*kc;
Gamma = (qa.*qc)./(qc.*(kj-ka)-qa.*(kj-kc));
Eta = qa./(qc-qa);