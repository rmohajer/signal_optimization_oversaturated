% plots for comparing the results for different Ns
fignum = 1;
loopnum = 1;

for flag_delay=[1,2,3]
    if flag_delay==1
        info_struct_QF = info_struct_QF_d;
        ds = 1; % to switch between the subplots
    elseif flag_delay==2
        info_struct_QF = info_struct_QF_th;
        ds=2;
    else
        info_struct_QF = info_struct_QF_sp;
        ds=3;
    end
    for N=Nset
        if info_struct_QF(N).exitflag >= 1 

            % back of queues evolution
            BQ1 = zeros(2*N,1);
            plottime1 = zeros(2*N,1);
            BQ2 = zeros(3*N,1);
            plottime2 = zeros(3*N,1);



            % green times
            GT1 = zeros(N,1);
            GT2 = zeros(N,1);

            for i=1:N
                BQ1((i-1)*2+1,1) = info_struct_QF(N).delta1(1,i);
                BQ1((i-1)*2+2,1) = info_struct_QF(N).delta3(1,i);
                plottime1((i-1)*2+1) = info_struct_QF(N).CritTime(1,(i-1)*4+1);
                plottime1((i-1)*2+2) = info_struct_QF(N).CritTime(1,(i-1)*4+3);


                BQ2((i-1)*3+1,1) = info_struct_QF(N).delta1(2,i);
                BQ2((i-1)*3+2,1) = info_struct_QF(N).delta2(2,i);
                BQ2((i-1)*3+3,1) = info_struct_QF(N).delta3(2,i);
                plottime2((i-1)*3+1) = info_struct_QF(N).CritTime(2,(i-1)*4+1);
                plottime2((i-1)*3+2) = info_struct_QF(N).CritTime(2,(i-1)*4+2);
                plottime2((i-1)*3+3) = info_struct_QF(N).CritTime(2,(i-1)*4+3);

                GT1(i,1) = info_struct_QF(N).Theta((i-1)*P+1);
                GT2(i,1) = info_struct_QF(N).Theta((i-1)*P+2);
            end



            % residual queues
            RQ1 = info_struct_QF(N).delta3(1,:);
            RQ2 = info_struct_QF(N).delta3(2,:);

            figure(fignum)
            subplot(2,3,ds)  
            plot(plottime1*60,BQ1*1000,plotStyle{loopnum})
            hold on
            xlabel('Critical time [min]')
            ylabel('Back of the queue at Approach 1[m]')
            legendInfo1{loopnum}= ['N = ' num2str(N)];
            legend(legendInfo1)

            subplot(2,3,ds+3)
            plot(plottime2*60,BQ2*1000,plotStyle{loopnum})
            hold on
            xlabel('Critical time [min]')
            ylabel('Back of the queue at Approach 2[m]')
            legend(legendInfo1)
            
            Delay = sum(info_struct_QF(N).delay,1)*60; % delay of each cycle in veh.min

            Mixed_OBJ = sum(info_struct_QF(N).Mixed_obj,1)*60; % mixed objective function
            
            Thr = sum(info_struct_QF(N).through,1); % throughput of each cycle in veh

            plottime = zeros(N,1);
            for i=1:N
                if i==1
                    plottime(i) = info_struct_QF(N).C(i);
                else

                    plottime(i) = plottime(i-1)+info_struct_QF(N).C(i);
                end
            end


            N_plot = N;

            Delay_CUM = cumsum(Delay);
            Mixed_OBJ_CUM = cumsum(Mixed_OBJ);
            Thr_CUM = cumsum(Thr);

            figure(fignum+1)
            subplot(3,3,ds)
            plot(plottime'*60,Delay_CUM,plotStyle{loopnum})
            hold on
            xlabel('Time [min]')
            ylabel('Accumulated total delay [veh.min]')
            legend(legendInfo1)

            subplot(3,3,ds+3)
            plot(plottime'*60,Mixed_OBJ_CUM,plotStyle{loopnum})
            hold on
            xlabel('Time [min]')
            ylabel('Accumulated mixed objective values')
            legend(legendInfo1)
            
            subplot(3,3,ds+6)
            plot(plottime'*60,Thr_CUM,plotStyle{loopnum})
            hold on
            xlabel('Time [min]')
            ylabel('Accumulated throughput [veh]')
            legend(legendInfo1)

            loopnum=loopnum+1;
       end
    end
    loopnum = 1;
end % for flag_delay



%% studying the variabiliy in delay and cycle time
Nplot = [];
delay_mean = [];
delay_std = [];
C_mean = [];
C_std = [];
for n=N1:N2
    
    if info_struct_QF_sp(n).exitflag>=1
        delay_n = sum(info_struct_QF_sp(n).delay)*60;
        info_struct_QF_sp(n).delay_mean = mean(delay_n);
        info_struct_QF_sp(n).delay_std = std(delay_n);
        Nplot=[Nplot,n];
        delay_mean = [delay_mean,info_struct_QF_sp(n).delay_mean];
        delay_std = [delay_std,info_struct_QF_sp(n).delay_std];
        
        C_mean = [C_mean, mean(info_struct_QF_sp(n).C)*60];
        C_std = [C_std, std(info_struct_QF_sp(n).C)*60];
    end % if exitflag>=1
    
    
end % for n

figure
subplot(1,2,1)
plot(Nplot,C_mean,'*-')
hold on
plot(Nplot,C_std,'ro-')
xlabel('Number of cycles to discharge the queues')
ylabel('Mean and Std. of cycle time per cycle [min]')
legend('mean cycle time', 'Std. of cycle time')

subplot(1,2,2)
plot(Nplot,delay_mean,'*-')
hold on
plot(Nplot,delay_std,'ro-')
xlabel('Number of cycles to discharge the queues')
ylabel('Mean and Std. of total vehicle delay per cycle [veh.min]')
legend('mean delay', 'Std. of delay')