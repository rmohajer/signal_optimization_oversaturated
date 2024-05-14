
% plots for comparing the results for different Ns
fignum = 1;
loopnum = 1;

for flag_delay=[1,0]
    if flag_delay
        info_struct = info_struct_d;
        ds = 1; % to switch between the subplots
    else
        info_struct = info_struct_th;
        ds=2;
    end
    
    
    
    for N=Nset
        if info_struct(N).exitflag == 1 

            % back of queues evolution
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

            if info_struct(N).exitflag_u==1
                for i=1:info_struct(N).Nu
                    BQ1(2*N+(i-1)*2+1,1) = info_struct(N).delta1_u(1,i);
                    BQ1(2*N+(i-1)*2+2,1) = info_struct(N).delta3_u(1,i);
                    plottime1(2*N+(i-1)*2+1) = info_struct(N).CritTime_u(1,(i-1)*4+1);
                    plottime1(2*N+(i-1)*2+2) = info_struct(N).CritTime_u(1,(i-1)*4+3);


                    BQ2(3*N+(i-1)*3+1,1) = info_struct(N).delta1_u(2,i);
                    BQ2(3*N+(i-1)*3+2,1) = info_struct(N).delta2_u(2,i);
                    BQ2(3*N+(i-1)*3+3,1) = info_struct(N).delta3_u(2,i);
                    plottime2(3*N+(i-1)*3+1) = info_struct(N).CritTime_u(2,(i-1)*4+1);
                    plottime2(3*N+(i-1)*3+2) = info_struct(N).CritTime_u(2,(i-1)*4+2);
                    plottime2(3*N+(i-1)*3+3) = info_struct(N).CritTime_u(2,(i-1)*4+3);


                    GT1(N+i,1) = info_struct(N).Theta_u((i-1)*P+1);
                    GT2(N+i,1) = info_struct(N).Theta_u((i-1)*P+2);
                end
            end

            % residual queues
            RQ1 = info_struct(N).delta3(1,:);
            RQ2 = info_struct(N).delta3(2,:);

            figure(fignum)
            subplot(2,2,ds)  
            plot(plottime1*60,BQ1*1000,plotStyle{loopnum})
            hold on
            xlabel('Critical time [min]')
            ylabel('Back of the queue at Approach 1[m]')
            legendInfo1{loopnum}= ['N = ' num2str(N)];
            legend(legendInfo1)

            subplot(2,2,ds+2)
            plot(plottime2*60,BQ2*1000,plotStyle{loopnum})
            hold on
            xlabel('Critical time [min]')
            ylabel('Back of the queue at Approach 2[m]')
            legend(legendInfo1)

            
            Delay = sum(info_struct(N).delay,1)*60; % delay of each cycle in veh.min
            Thr = sum(info_struct(N).through,1); % throughput of each cycle in veh
            plottime = zeros(N,1);
            for i=1:N
                if i==1
                    plottime(i) = info_struct(N).C(i);
                else

                    plottime(i) = plottime(i-1)+info_struct(N).C(i);
                end
            end
            N_plot = N;
            if info_struct(N).exitflag_u==1
                Delay = [Delay, sum(info_struct(N).delay_u,1)*60];
                Thr = [Thr, sum(info_struct(N).through_u,1)];
                N_plot = N+info_struct(N).Nu;

                for i=1:info_struct(N).Nu
                    plottime(N+i) = plottime(N+i-1)+info_struct(N).C_u(i);

                end
            end

            Delay_CUM = cumsum(Delay);
            Thr_CUM = cumsum(Thr);

            figure(fignum+1)
            
            subplot(2,2,ds)       
            plot(plottime'*60,Delay_CUM,plotStyle{loopnum})
            hold on
            xlabel('Time [min]')
            ylabel('Accumulated total delay [veh.min]')
            legend(legendInfo1)
            
            subplot(2,2,ds+2)            
            plot(plottime'*60,Thr_CUM,plotStyle{loopnum})
            hold on
            xlabel('Time [min]')
            ylabel('Accumulated Throughput [veh]')
            legend(legendInfo1)



            loopnum=loopnum+1;
       end
    end
    loopnum = 1;
end % for flag_delay
    

%% studying the variabiliy in delay and cycle time

figure
for flag_delay=[1,0]
    if flag_delay
        info_struct = info_struct_d;
        ds = 0; % to switch between the subplots
    else
        info_struct = info_struct_th;
        ds=3;
    end

    Nplot = [];
    delay_mean = [];
    delay_std = [];
    duration = [];
    C_mean = [];
    C_std = [];
    for n=N1:N2

        if info_struct(n).exitflag>=1
            delay_n = sum(info_struct(n).delay)*60;
            info_struct(n).delay_mean = mean(delay_n);
            info_struct(n).delay_std = std(delay_n);
            Nplot=[Nplot,n];
            delay_mean = [delay_mean,info_struct(n).delay_mean];
            delay_std = [delay_std,info_struct(n).delay_std];
            duration = [duration, info_struct(n).solvedur*1000];

            C_mean = [C_mean, mean(info_struct(n).C)*60];
            C_std = [C_std, std(info_struct(n).C)*60];
        end % if exitflag>=1


    end % for n

    
    subplot(2,3,ds+1)
    plot(Nplot,C_mean,'*-')
    hold on
    plot(Nplot,C_std,'ro-')
    xlabel('Number of cycles to discharge the queues')
    ylabel('Mean and Std. of cycle time per cycle [min]')
    legend('mean cycle time', 'Std. of cycle time')

    subplot(2,3,ds+2)
    plot(Nplot,delay_mean,'*-')
    hold on
    plot(Nplot,delay_std,'ro-')
    xlabel('Number of cycles to discharge the queues')
    ylabel('Mean and Std. of total vehicle delay per cycle [veh.min]')
    legend('mean delay', 'Std. of delay')
    
    subplot(2,3,ds+3)
    plot(Nplot,duration,'*-')
    hold on
    xlabel('Number of cycles to discharge the queues')
    ylabel('Solution time [ms]')
end % for flag_delay




%% plotting the queuing theory
% loopnum = 1;
% fignum = fignum+3;
% for N=Nset
%     
%    
%     if info_struct_q(N).exitflag == 1 
% 
%         % back of the queues based on the queuing theory model
%         BQq1 = zeros(2*N,1);
%         plottimeq1 = zeros(2*N,1);
%         BQq2 = zeros(3*N,1);
%         plottimeq2 = zeros(3*N,1);
%         
%         % back of the queues based on the LWR theory model
%         BQqLWR1 = zeros(2*N,1);
%         BQqLWR2 = zeros(2*N,1);
%         
%         % green times
%         GT1_q = zeros(N,1);
%         GT2_q = zeros(N,1);
% 
%         for i=1:N
%             BQq1((i-1)*2+1,1) = info_struct_q(N).delta1_q(1,i)/kj(1);
%             BQq1((i-1)*2+2,1) = info_struct_q(N).delta3_q(1,i)/kj(1);
%             plottimeq1((i-1)*2+1) = info_struct_q(N).CritTime(1,(i-1)*4+1);
%             plottimeq1((i-1)*2+2) = info_struct_q(N).CritTime(1,(i-1)*4+3);
% 
%             BQqLWR1((i-1)*2+1,1) = info_struct_q(N).delta1(1,i);
%             BQqLWR1((i-1)*2+2,1) = info_struct_q(N).delta3(1,i);
%             
%             
%             BQq2((i-1)*3+1,1) = info_struct_q(N).delta1_q(2,i)/kj(2);
%             BQq2((i-1)*3+2,1) = info_struct_q(N).delta2_q(2,i)/kj(2);
%             BQq2((i-1)*3+3,1) = info_struct_q(N).delta3_q(2,i)/kj(2);
%             plottimeq2((i-1)*3+1) = info_struct_q(N).CritTime(2,(i-1)*4+1);
%             plottimeq2((i-1)*3+2) = info_struct_q(N).CritTime(2,(i-1)*4+2);
%             plottimeq2((i-1)*3+3) = info_struct_q(N).CritTime(2,(i-1)*4+3);
%             
%             
%             BQqLWR2((i-1)*3+1,1) = info_struct_q(N).delta1(2,i);
%             BQqLWR2((i-1)*3+2,1) = info_struct_q(N).delta2(2,i);
%             BQqLWR2((i-1)*3+3,1) = info_struct_q(N).delta3(2,i);
% 
%             GT1_q(i,1) = info_struct_q(N).Theta((i-1)*P+1);
%             GT2_q(i,1) = info_struct_q(N).Theta((i-1)*P+2);
%         end
% 
% 
%         % residual queues
%         RQq1 = info_struct_q(N).delta3_q(1,:)/kj(1);
%         RQq2 = info_struct_q(N).delta3_q(2,:)/kj(2);
%         
%         % residual queues based on the LWR theory model
%         RQqLWR1 = info_struct_q(N).delta3(1,:);
%         RQqLWR2 = info_struct_q(N).delta3(2,:);
%         
% 
%         Delayq = sum(info_struct_q(N).delay_q,1)*60; % delay of each cycle in veh.min
%         DelayqLWR = sum(info_struct_q(N).delay,1)*60; % delay of each cycle in veh.min
%         
%         plottimeq = zeros(N,1);
%         for i=1:N
%             if i==1
%                 plottimeq(i) = info_struct_q(N).C(i);
%             else
% 
%                 plottimeq(i) = plottimeq(i-1)+info_struct_q(N).C(i);
%             end
%         end
%         Delay_CUMq = cumsum(Delayq);
%         Delay_CUMqLWR = cumsum(DelayqLWR);
% 
%         figure(fignum)
%         subplot(4,2,1)  
%         plot(plottimeq1*60,BQq1*1000,plotStyle{loopnum})
%         hold on
%         xlabel('Critical time [min]')
%         ylabel('Back of the queue at Approach 1[m] (QTB model)')
%         legendInfo1{loopnum}= ['N = ' num2str(N)];
%         legend(legendInfo1)
%         
%         subplot(4,2,2)  
%         plot(plottimeq1*60,BQqLWR1*1000,plotStyle_com{loopnum})
%         hold on
%         xlabel('Critical time [min]')
%         ylabel('Back of the queue at Approach 1[m] (LWRTB model)')
%         legendInfo1{loopnum}= ['N = ' num2str(N)];
%         legend(legendInfo1)
% 
%         subplot(4,2,3)
%         plot(plottimeq2*60,BQq2*1000,plotStyle{loopnum})
%         hold on
%         xlabel('Critical time [min]')
%         ylabel('Back of the queue at Approach 2 [m] (QTB model)')
%         legend(legendInfo1)
%         
%         subplot(4,2,4)
%         plot(plottimeq2*60,BQqLWR2*1000,plotStyle_com{loopnum})
%         hold on
%         xlabel('Critical time [min]')
%         ylabel('Back of the queue at Approach 2 [m] (LWRTB model)')
%         legend(legendInfo1)
% 
%         subplot(4,2,5)
%         plot(1:N, RQq1*1000,plotStyle{loopnum})
%         hold on
%         xlabel('Cycle number')
%         ylabel('Residual Queue at Approach 1[m] (QTB model)')
%         legendInfo1{loopnum}= ['N = ' num2str(N)];
%         
%         subplot(4,2,6)
%         plot(1:N, RQqLWR1*1000,plotStyle_com{loopnum})
%         hold on
%         xlabel('Cycle number')
%         ylabel('Residual Queue at Approach 1[m] (LWRTB model)')
%         legendInfo1{loopnum}= ['N = ' num2str(N)];
% 
%         subplot(4,2,7)
%         plot(1:N, RQq2*1000,plotStyle{loopnum})
%         hold on
%         xlabel('Cycle number')
%         ylabel('Residual Queue at Approach 2[m] (QTB model)')
%         legendInfo1{loopnum}= ['N = ' num2str(N)];
%         
%         subplot(4,2,8)
%         plot(1:N, RQqLWR2*1000,plotStyle_com{loopnum})
%         hold on
%         xlabel('Cycle number')
%         ylabel('Residual Queue at Approach 2[m] (LWRTB model)')
%         legendInfo1{loopnum}= ['N = ' num2str(N)];
% 
% 
%         figure(fignum+1)
%         subplot(2,2,1)
%         plot(1:N,Delayq,plotStyle{loopnum})
%         hold on
%         xlabel('Cycle number')
%         ylabel('Total delay [veh.min] (QTB model)')
%         legend(legendInfo1)
%         
%         
%         subplot(2,2,2)
%         plot(1:N,DelayqLWR,plotStyle_com{loopnum})
%         hold on
%         xlabel('Cycle number')
%         ylabel('Total delay [veh.min] (LWRTB model)')
%         legend(legendInfo1)
% 
%         subplot(2,2,3)
%         plot(plottimeq'*60,Delay_CUMq,plotStyle{loopnum})
%         hold on
%         xlabel('Time [min]')
%         ylabel('Accumulated total delay [veh.min] (QTB model)')
%         legend(legendInfo1)
%         
%         subplot(2,2,4)
%         plot(plottimeq'*60,Delay_CUMqLWR,plotStyle_com{loopnum})
%         hold on
%         xlabel('Time [min]')
%         ylabel('Accumulated total delay [veh.min] (LWRTB model)')
%         legend(legendInfo1)
%         
%         loopnum=loopnum+1;
%     end % info_struct_q(N).exitflag == 1 
%         
% end


%% Plotting the results for 1 cycle discharge

% eta1 = info_struct(1).Eta(1,1);
% eta2 = info_struct(1).Eta(2,1);
% g1 = 0:g1max/10:g1max;
% 
% g2_qd = eta2*g1+eta2*l(1)+kj(2)/qc(2)*delt_0(2);
% 
% g2_maxc = -g1+cmax-L;
% 
% g2_1 = eta2*g1+eta2*L;
% g2_2 = 1/eta1*g1-L;
% 
% Nbb = sum(info_structBB.Nopt);
% 
% N=Nbb;
% Thx_N = zeros(1,N);
% Thy_N = zeros(1,N);
% for i=1:N
%     Thx_N(i) = info_struct(N).Theta(2*i-1);
%     Thy_N(i) = info_struct(N).Theta(2*i);
% end
% 
% eta1_N = min(info_struct(N).Eta(1,:));
% eta2_N = min(info_struct(N).Eta(2,:));
% 
% g2_N1 = eta2_N*g1+eta2_N*L; % for undersaturation constraint
% g2_N2 = 1/eta1_N*g1-L;% for undersaturation constraint
% 
% 
% Thxbb = zeros(1,Nbb);
% Thybb = zeros(1,Nbb);
% for i=1:Nbb
%     Thxbb(i) = info_structBB.Theta(2*i-1);
%     Thybb(i) = info_structBB.Theta(2*i);
% end
% 
% 
% eta1_bb = min(info_structBB.Eta(1,:));
% eta2_bb = min(info_structBB.Eta(2,:));
% 
% g2_bb1 = eta2_bb*g1+eta2_bb*L; % for undersaturation constraint
% g2_bb2 = 1/eta1_bb*g1-L;% for undersaturation constraint
% 
% figure
% 
% plot([eta1*L*3600,g1max*3600],[0,(1/eta1*g1max-L)*3600],'b')
% hold on
% plot(g1*3600,g2_1*3600,'b')
% plot(info_struct(1).Theta(1)*3600,info_struct(1).Theta(2)*3600,'*')
% plot(g1*3600,g2_qd*3600,'r--')
% plot([kj(1)/qc(1)*delt_0(1)*3600,kj(1)/qc(1)*delt_0(1)*3600],[0,g1max*3600],'r--')
% 
% plot([eta1_N*L*3600,g1max*3600],[0,(1/eta1_N*g1max-L)*3600],'m')
% plot(g1*3600,g2_N1*3600,'m')
% plot(Thx_N*3600,Thy_N*3600,'*m')
% 
% plot([eta1_bb*L*3600,g1max*3600],[0,(1/eta1_bb*g1max-L)*3600],'k')
% plot(g1*3600,g2_bb1*3600,'k')
% plot(Thxbb*3600,Thybb*3600,'*k')
% 
% plot(g1*3600,g2_maxc*3600,'g--')    


