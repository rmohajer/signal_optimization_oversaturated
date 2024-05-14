randq1 = randi([qmin(1),qmax(1)],1,T_total/ts); % random flows of approach 1
randq2 = randi([qmin(2),qmax(2)],1,T_total/ts); % random flows of approach 2

for i=1:size(randq1,1)
    if randq1(i)>1000 && randq2(i)>730
        randq1(i) = 1000;
        randq2(i) = 730;
    elseif randq1(i)/qc1+randq2(i)/qc2+L/cmax>1
        if randq1(i)>1000
            randq1(i)=1000;
        else
            randq2(i) = 730;
        end
    end
end


for i=2:size(randq1,2)
    %if there is a big jump in the flow, set it to a threshold of 300 veh/h/l
    % Approach 1
   if randq1(i)-randq1(i-1)>200
       randq1(i)=randq1(i-1)+100;
   elseif randq1(i)-randq1(i-1)<-200
        randq1(i)=randq1(i-1)-100;
   end
   % Approach 2
   if randq2(i)-randq2(i-1)>200
       randq2(i)=randq2(i-1)+100;
   elseif randq2(i)-randq2(i-1)<-200
        randq2(i)=randq2(i-1)-100;
   end
end

Qa = [[randq1,randq1(end)];[randq2,randq2(end)]] %veh/h


T_cycle = 0:ts:T_total; % min

figure(1)
stairs((0:ts:T_total)*60,Qa(1,:))
hold on
stairs((0:ts:T_total)*60,Qa(2,:),'r')
xlabel('Time [min]')
ylabel('q^a [veh/h]')
legend('Approach 1','Approach 2')