cyc_num = T_total/ts;
oversat = zeros(1,cyc_num);

while prod(oversat)<1
    randq1 = randi([qmin(1),qmax(1)],1,T_generate/ts); % random flows of approach 1
    randq2 = randi([qmin(2),qmax(2)],1,T_generate/ts); % random flows of approach 2



    for i=2:cyc_num
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

    oversat = zeros(1,cyc_num);
    for cyc=1:cyc_num
        oversat(cyc) = Qa(1,cyc)/qc(1)+Qa(2,cyc)/qc(2)+L/cmax>1;
    end
end

disp('Cycles that are oversaturated')
oversat

T_cycle = 0:ts:T_generate; % min

figure(1)
stairs((0:ts:T_generate)*60,Qa(1,:))
hold on
stairs((0:ts:T_generate)*60,Qa(2,:),'r')
xlabel('Time [min]')
ylabel('q^a [veh/h]')
legend('Approach 1','Approach 2')