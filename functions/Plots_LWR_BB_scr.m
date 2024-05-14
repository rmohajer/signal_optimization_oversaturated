% plots for comparing the bang bang results

Nset = [NBB-1,NBB,NBB+1];
loopnum = 1;
fignum = fignum+3;

if NBB<Nmax
    for N=Nset

        %first plot the bang bang results
        if N==NBB 

            % Back of the queues
            BQBB1 = zeros(2*N,1);
            plottimeBB1 = zeros(2*N,1);
            BQBB2 = zeros(3*N,1);
            plottimeBB2 = zeros(3*N,1);
            
            % green times
            GTBB1 = zeros(N,1);
            GTBB2 = zeros(N,1);

            for i=1:N
                BQBB1((i-1)*2+1,1) = info_structBB.delta1(1,i);
                BQBB1((i-1)*2+2,1) = info_structBB.delta3(1,i);
                
                plottimeBB1((i-1)*2+1) = info_structBB.CritTime(1,(i-1)*4+1);
                plottimeBB1((i-1)*2+2) = info_structBB.CritTime(1,(i-1)*4+3);
               
                BQBB2((i-1)*3+1,1) = info_structBB.delta1(2,i);
                BQBB2((i-1)*3+2,1) = info_structBB.delta2(2,i);
                BQBB2((i-1)*3+3,1) = info_structBB.delta3(2,i);
                plottimeBB2((i-1)*3+1) = info_structBB.CritTime(2,(i-1)*4+1);
                plottimeBB2((i-1)*3+2) = info_structBB.CritTime(2,(i-1)*4+2);
                plottimeBB2((i-1)*3+3) = info_structBB.CritTime(2,(i-1)*4+3);
                
                GTBB1(i,1) = info_structBB.Theta((i-1)*P+1);
                GTBB2(i,1) = info_structBB.Theta((i-1)*P+2);
            end
            
            
            % residual queues
            RQBB1 = info_structBB.delta3(1,:);
            RQBB2 = info_structBB.delta3(2,:);

            figure(fignum)
            subplot(2,2,1)  
            plot(plottimeBB1*60,BQBB1*1000,'k*-')
            hold on
            xlabel('Critical time [min]')
            ylabel('Back of the queue at Approach 1[m]')
            legendInfo1{loopnum}= ['N = ' num2str(N),' (bang bang)'];
            legend(legendInfo1)

            subplot(2,2,2)
            plot(plottimeBB2*60,BQBB2*1000,'k*-')
            hold on
            xlabel('Critical time [min]')
            ylabel('Back of the queue at Approach 2[m]')
            legend(legendInfo1)

            subplot(2,2,3)
            plot(1:N, RQBB1*1000,plotStyle{loopnum})
            hold on
            xlabel('Cycle number')
            ylabel('Residual Queue at Approach 1[m]')
            legendInfo1{loopnum}= ['N = ' num2str(N),' (bang bang)'];

            subplot(2,2,4)
            plot(1:N, RQBB2*1000,'k*-')
            hold on
            xlabel('Cycle number')
            ylabel('Residual Queue at Approach 2[m]')
            legendInfo1{loopnum}= ['N = ' num2str(N),' (bang bang)'];

            Delay = sum(info_structBB.delay,1)*60; % delay of each cycle in veh.min
            plottime = zeros(N,1);
            for i=1:N
                if i==1
                    plottime(i) = info_structBB.C(i);
                else

                    plottime(i) = plottime(i-1)+info_structBB.C(i);
                end
            end

            Delay_CUM = cumsum(Delay);

            figure(fignum+1)
            subplot(1,2,1)
            plot(1:N,Delay,'k*-')
            hold on
            xlabel('Cycle number')
            ylabel('Total delay [veh.min]')
            legend(legendInfo1)

            subplot(1,2,2)
            plot(plottime'*60,Delay_CUM,'k*-')
            hold on
            xlabel('Time [min]')
            ylabel('Accumulated total delay [veh.min]')
            legend(legendInfo1)
            loopnum=loopnum+1;
        end

        % now plot the normal algorithm results

        if info_struct(N).exitflag == 1 

            % residual queues
            BQ1 = zeros(2*N,1);
            plottime1 = zeros(2*N,1);
            BQ2 = zeros(3*N,1);
            plottime2 = zeros(3*N,1);
            
            % green times
            GT1 = zeros(N,1);
            GT2 = zeros(N,1);

            for i=1:N
                BQ1((i-1)*2+1,1) = info_struct(N).delta1(1,i);
                BQ1((i-1)*2+2,1) = info_struct(N).delta3(1,i);
                plottime1((i-1)*2+1) = info_struct(N).CritTime(1,(i-1)*4+1);
                plottime1((i-1)*2+2) = info_struct(N).CritTime(1,(i-1)*4+3);


                BQ2((i-1)*3+1,1) = info_struct(N).delta1(2,i);
                BQ2((i-1)*3+2,1) = info_struct(N).delta2(2,i);
                BQ2((i-1)*3+3,1) = info_struct(N).delta3(2,i);
                plottime2((i-1)*3+1) = info_struct(N).CritTime(2,(i-1)*4+1);
                plottime2((i-1)*3+2) = info_struct(N).CritTime(2,(i-1)*4+2);
                plottime2((i-1)*3+3) = info_struct(N).CritTime(2,(i-1)*4+3);
                
                GT1(i,1) = info_struct(N).Theta((i-1)*P+1);
                GT2(i,1) = info_struct(N).Theta((i-1)*P+2);
            end

%             if info_struct(N).exitflag_u==1
%                 for i=1:info_struct(N).Nu
%                     BQ1(2*N+(i-1)*2+1,1) = info_struct(N).delta1_u(1,i);
%                     BQ1(2*N+(i-1)*2+2,1) = info_struct(N).delta3_u(1,i);
%                     plottime1(2*N+(i-1)*2+1) = info_struct(N).CritTime_u(1,(i-1)*4+1);
%                     plottime1(2*N+(i-1)*2+2) = info_struct(N).CritTime_u(1,(i-1)*4+3);
% 
% 
%                     BQ2(3*N+(i-1)*3+1,1) = info_struct(N).delta1_u(2,i);
%                     BQ2(3*N+(i-1)*3+2,1) = info_struct(N).delta2_u(2,i);
%                     BQ2(3*N+(i-1)*3+3,1) = info_struct(N).delta3_u(2,i);
%                     plottime2(3*N+(i-1)*3+1) = info_struct(N).CritTime_u(2,(i-1)*4+1);
%                     plottime2(3*N+(i-1)*3+2) = info_struct(N).CritTime_u(2,(i-1)*4+2);
%                     plottime2(3*N+(i-1)*3+3) = info_struct(N).CritTime_u(2,(i-1)*4+3);
%                     
%                     GT1(N+i,1) = info_struct(N).Theta_u((i-1)*P+1);
%                     GT2(N+i,1) = info_struct(N).Theta_u((i-1)*P+2);
%                 end
%             end
            
            % residual queues
            RQ1 = info_struct(N).delta3(1,:);
            RQ2 = info_struct(N).delta3(2,:);
            
            figure(fignum)
            subplot(2,2,1)  
            plot(plottime1*60,BQ1*1000,plotStyle{loopnum})
            hold on
            xlabel('Critical time [min]')
            ylabel('Back of the queue at Approach 1[m]')
            legendInfo1{loopnum}= ['N = ' num2str(N)];
            legend(legendInfo1)

            subplot(2,2,2)
            plot(plottime2*60,BQ2*1000,plotStyle{loopnum})
            hold on
            xlabel('Critical time [min]')
            ylabel('Back of the queue at Approach 2[m]')
            legend(legendInfo1)

            subplot(2,2,3)
            plot(1:N, RQ1*1000,plotStyle{loopnum})
            hold on
            xlabel('Cycle number')
            ylabel('Residual Queue at Approach 1[m]')
            legendInfo1{loopnum}= ['N = ' num2str(N)];

            subplot(2,2,4)
            plot(1:N, RQ2*1000,plotStyle{loopnum})
            hold on
            xlabel('Cycle number')
            ylabel('Residual Queue at Approach 2[m]')
            legendInfo1{loopnum}= ['N = ' num2str(N)];


            Delay = sum(info_struct(N).delay,1)*60; % delay of each cycle in veh.min
            plottime = zeros(N,1);
            for i=1:N
                if i==1
                    plottime(i) = info_struct(N).C(i);
                else

                    plottime(i) = plottime(i-1)+info_struct(N).C(i);
                end
            end

            Delay_CUM = cumsum(Delay);

            figure(fignum+1)
            subplot(1,2,1)
            plot(1:N,Delay,plotStyle{loopnum})
            hold on
            xlabel('Cycle number')
            ylabel('Total delay [veh.min]')
            legend(legendInfo1)

            subplot(1,2,2)
            plot(plottime'*60,Delay_CUM,plotStyle{loopnum})
            hold on
            xlabel('Time [min]')
            ylabel('Accumulated total delay [veh.min]')
            legend(legendInfo1)

            loopnum=loopnum+1;
       end 

    end
    
else
    
   % residual queues
    N=NBB;
   
    RQBB1 = zeros(2*N,1);
    plottimeBB1 = zeros(2*N,1);
    RQBB2 = zeros(3*N,1);
    plottimeBB2 = zeros(3*N,1);

    for i=1:N
        RQBB1((i-1)*2+1,1) = info_structBB.delta1(1,i);
        RQBB1((i-1)*2+2,1) = info_structBB.delta3(1,i);
        plottimeBB1((i-1)*2+1) = info_structBB.CritTime(1,(i-1)*4+1);
        plottimeBB1((i-1)*2+2) = info_structBB.CritTime(1,(i-1)*4+3);


        RQBB2((i-1)*3+1,1) = info_structBB.delta1(2,i);
        RQBB2((i-1)*3+2,1) = info_structBB.delta2(2,i);
        RQBB2((i-1)*3+3,1) = info_structBB.delta3(2,i);
        plottimeBB2((i-1)*3+1) = info_structBB.CritTime(2,(i-1)*4+1);
        plottimeBB2((i-1)*3+2) = info_structBB.CritTime(2,(i-1)*4+2);
        plottimeBB2((i-1)*3+3) = info_structBB.CritTime(2,(i-1)*4+3);
    end

    figure(fignum)
    subplot(1,2,1)  
    plot(plottimeBB1*60,RQBB1*1000,plotStyle{loopnum})
    hold on


    xlabel('Critical time [min]')
    ylabel('Residual Queue at Approach 1[m]')
    legendInfo1{loopnum}= ['N = ' num2str(N),' (bang bang)'];
    legend(legendInfo1)

    subplot(1,2,2)
    plot(plottimeBB2*60,RQBB2*1000,plotStyle{loopnum})
    hold on
    xlabel('Critical time [min]')
    ylabel('Residual Queue at Approach 2[m]')
    legend(legendInfo1)


    Delay = sum(info_structBB.delay,1)*60; % delay of each cycle in veh.min
    plottime = zeros(N,1);
    for i=1:N
        if i==1
            plottime(i) = info_structBB.C(i);
        else

            plottime(i) = plottime(i-1)+info_structBB.C(i);
        end
    end

    Delay_CUM = cumsum(Delay);

    figure(fignum+1)
    subplot(1,2,1)
    plot(1:N,Delay,plotStyle{loopnum})
    hold on
    xlabel('Cycle number')
    ylabel('Total delay [veh.min]')
    legend(legendInfo1)

    subplot(1,2,2)
    plot(plottime'*60,Delay_CUM,plotStyle{loopnum})
    hold on
    xlabel('Time [min]')
    ylabel('Accumulated total delay [veh.min]')
    legend(legendInfo1)
    loopnum=loopnum+1; 
    
end