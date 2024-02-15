%dl  = [1, 0.05, 0.005];
%delta = 0.95;
test = [7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140];
net = ["Twaren.graphml","Abvt.graphml","GtsHungary.graphml","GtsSlovakia.graphml","Intranetwork.graphml","BtNorthAmerica.graphml","Bellsouth.graphml","Iris.graphml","Renater2010.graphml","GtsCe.graphml"];
%test = [7,14];
%test = 42;
prob_success_mmtc = zeros(10,4);
prob_success_embb = zeros(10,4);
prob_success_urllc = zeros(10,4);
power_opt = zeros(10,length(net));
power_sta = zeros(10,length(net));
power_prop_50 = zeros(10,length(net));
power_prop_75 = zeros(10,length(net));
power_prop_100= zeros(10,length(net));
power_opt_final = zeros(10,1);
num_slices = 3;
power_slice_final = zeros(10,num_slices);
err = zeros(10,1);
power_sta_final = zeros(10,1);
power_prop_50_final = zeros(10,1);
power_prop_75_final = zeros(10,1);
power_prop_100_final= zeros(10,1);
%arrival_set = [1,10,20,30,40,50,60,70,80,90,100,200,300,400,500];
%com_mmtc = [];
%com_embb = [];
%com_urllc = [];
%prob_fail = [0.1,0.01,0.007];
scale_factor = [0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2];
%for iter = 1
    %i = arrival_set(iter);
    %flag = [100 20 5];
    %for s = 1:num_slices
    %    arrival_rate(s) = i*flag(s);
    %end
    %arrival_mmtc(iter) = flag(1)*i;
    %arrival_embb(iter) = flag(2)*i;
    %arrival_urllc(iter) = flag(3)*i;
   % active_set = [1,2,3];
    %arrival_set = [arrival_set i*10];
    prob_50 = {};
    prob_75 = {};
    prob_100 = {};
    prob_opt ={};
    prob_opt_matrix = {};
    prob_prop_50_new = {};
    prob_prop_75_new = {};
    prob_prop_100_new = {};
    prob_prop_50 = {};
    prob_prop_75 = {};
    prob_prop_100 = {};
    count_all = 0;
    count_drop = 0;
    iter = 5;
    for net_number = 1:10
    for seed_chose = 1:length(test)
    power_opt_stochastic = [];
    prob_opt_stochastic = [];
    seed = test(seed_chose);
    for iteration = 1:5
    filename = ['final_placement_1s_50ms_5ms_slices_final_',convertStringsToChars(net(net_number)),'_',num2str(seed),'_complexity_delay1e4_iter200_contolerance1-2_inter_5_test_greedy_',num2str(iteration),'.mat'];
    load(filename);
    count_all = count_all + 1;
    power_opt_stochastic(iteration) = sum(f_opt{seed}(iter,:));
    prob_opt_stochastic(iteration,:) = prob_opt{seed}(iter,:);
        if sum(isnan(prob_opt_stochastic(iteration,:))) ~= 0
            power_opt_stochastic(iteration) = 1e6;
            count_drop = count_drop + 1;
        end
        if sum(isnan(prob_prop_50{seed}(iter,:))) == 0
            prob_prop_50_new{seed}(iter,:) = prob_prop_50{seed}(iter,:); 
        end
        if sum(isnan(prob_prop_75{seed}(iter,:))) == 0
            prob_prop_75_new{seed}(iter,:) = prob_prop_75{seed}(iter,:); 
        end
        if sum(isnan(prob_prop_100{seed}(iter,:))) == 0
            prob_prop_100_new{seed}(iter,:) = prob_prop_100{seed}(iter,:); 
        end
    end
    [val,pos] = min(power_opt_stochastic);
    filename = ['final_placement_1s_50ms_5ms_slices_final_',convertStringsToChars(net(net_number)),'_',num2str(seed),'_complexity_delay1e4_iter200_contolerance1-2_inter_5_test_greedy_',num2str(pos),'.mat'];
    filename
    load(filename);
    %%% propotional strategy:
    %[f_prop_50_new{seed}(iter,:),f_prop_75_new{seed}(iter,:),f_prop_100_new{seed}(iter,:),prob_prop_50_new{seed}(iter,:),prob_prop_75_new{seed}(iter,:),prob_prop_100_new{seed}(iter,:)] = load_propstra_com_final(seed,iter);
    %f_prop_50{7}(1,:)
    %prob_prop_50{seed}(iter,:) = prob_prop_50_new{seed}(iter,:);
    %prob_prop_75{seed}(iter,:) = prob_prop_75_new{seed}(iter,:);
    %prob_prop_100{seed}(iter,:) = prob_prop_100_new{seed}(iter,:);
    t{iter}(seed_chose,:) = t_app;
    prob_opt_matrix{seed} = (prob_opt{seed}>=0) & (prob_opt{seed}<=1.0001);
    prob_success_mmtc(net_number,1) = prob_success_mmtc(net_number,1) + prob_opt{seed}(iter,1)*prob_opt_matrix{seed}(iter,1);
    prob_success_embb(net_number,1) = prob_success_embb(net_number,1) + prob_opt{seed}(iter,2)*prob_opt_matrix{seed}(iter,2);
    prob_success_urllc(net_number,1) = prob_success_urllc(net_number,1) + prob_opt{seed}(iter,3)*prob_opt_matrix{seed}(iter,3);
    
    prob_50{seed} = (prob_prop_50_new{seed}>=0) & (prob_prop_50_new{seed}<=1.0001);
    if prob_50{seed}(iter,1) == 1
    prob_success_mmtc(net_number,2) = prob_success_mmtc(net_number,2) + prob_prop_50_new{seed}(iter,1)*prob_50{seed}(iter,1);
    end
    if prob_50{seed}(iter,2) == 1
    prob_success_embb(net_number,2) = prob_success_embb(net_number,2) + prob_prop_50_new{seed}(iter,2)*prob_50{seed}(iter,2);
    end
    if prob_50{seed}(iter,3) == 1
    prob_success_urllc(net_number,2) = prob_success_urllc(net_number,2) + prob_prop_50_new{seed}(iter,3)*prob_50{seed}(iter,3);
    end
    prob_75{seed} = (prob_prop_75_new{seed}>=0) & (prob_prop_75_new{seed}<=1.0001);
    if prob_75{seed}(iter,1) == 1
    prob_success_mmtc(net_number,3) = prob_success_mmtc(net_number,3) + prob_prop_75_new{seed}(iter,1)*prob_75{seed}(iter,1);
    end
    if prob_75{seed}(iter,2) == 1
    prob_success_embb(net_number,3) = prob_success_embb(net_number,3) + prob_prop_75_new{seed}(iter,2)*prob_75{seed}(iter,2);
    end
    if prob_75{seed}(iter,3) == 1
    prob_success_urllc(net_number,3) = prob_success_urllc(net_number,3) + prob_prop_75_new{seed}(iter,3)*prob_75{seed}(iter,3);
    end
    
    prob_100{seed} = (prob_prop_100_new{seed}>=0) & (prob_prop_100_new{seed}<=1.0001);
    if prob_100{seed}(iter,1) == 1
    prob_success_mmtc(net_number,4) = prob_success_mmtc(net_number,4) + prob_prop_100_new{seed}(iter,1)*prob_100{seed}(iter,1);
    end
    if prob_100{seed}(iter,2) == 1
    prob_success_embb(net_number,4) = prob_success_embb(net_number,4) + prob_prop_100_new{seed}(iter,2)*prob_100{seed}(iter,2);
    end
    if prob_100{seed}(iter,3) == 1
    prob_success_urllc(net_number,4) = prob_success_urllc(net_number,4) + prob_prop_100_new{seed}(iter,3)*prob_100{seed}(iter,3);
    end
    power_opt(net_number,seed_chose) = sum(f_opt{seed}(iter,:));
    power_sta(net_number,seed_chose) = sum(f_sta{seed}(iter,:));
    power_prop_50(net_number,seed_chose) = sum(f_prop_50{seed}(iter,:));  
    power_prop_75(net_number,seed_chose) = sum(f_prop_75{seed}(iter,:));  
    power_prop_100(net_number,seed_chose) = sum(f_prop_100{seed}(iter,:));  
        for n = 1:num_slices
            power_slice_final(net_number,n) = power_slice_final(net_number,n) + f_opt{seed}(iter, n);
        end
    end
    prob_success_mmtc(net_number,:) = prob_success_mmtc(net_number,:)/length(test);
    prob_success_embb(net_number,:) = prob_success_embb(net_number,:)/length(test);
    prob_success_urllc(net_number,:) = prob_success_urllc(net_number,:)/length(test);
    power_opt_final(net_number) = mean(power_opt(net_number,:));
    power_sta_final(net_number) = mean(power_sta(net_number,:));
    power_prop_50_final(net_number) = mean(power_prop_50(net_number,:));
    power_prop_75_final(net_number) = mean(power_prop_75(net_number,:));
    power_prop_100_final(net_number) = mean(power_prop_100(net_number,:));
    [muhat,sigmahat,muci,sigmaci] = normfit(power_opt(net_number,:));
    errneg(net_number) = muci(1);
    errpos(net_number) = muci(2);

    end
    power_slice_final =   power_slice_final /length(test);
fprintf( 'drop percentage: %f \n', count_drop/count_all);
y = [prob_success_mmtc(1,:);prob_success_embb(1,:);prob_success_urllc(1,:)];
%y2 = [power_sta,power_opt,power_prop_50,power_prop_100,]';
figure;
%h1 = bar(scale_factor,y(1));
%ind = [1:length(arrival_set)];
%h1 = bar(scale_factor,prob_success_mmtc);
SH = shadowHist(prob_success_mmtc,'ShadowType',{'/','.','x','|'});
SH = SH.draw;
SH = SH.legend({'OptRes','PropRes50','PropRes75','PropRes100'});
%xlabel('Arrival rate of slice mMTC (\cdot Default value) [Reqs/sec]');
xlabel('Duplication probability of slice mMTC (\cdot Default value) [instr/req, bytes/req]');
%xlabel('Routing probability of slice mMTC');
ylabel('Reliability (mMTC)');
%legend('OptRes','PropRes50','PropRes75','PropRes100','Location','best');
%set(gca,'XTickLabel',{'mMTC','eMBB','URLLC'});
figure;
%h2 = bar(scale_factor,prob_success_embb);
SH = shadowHist(prob_success_embb,'ShadowType',{'/','.','x','|'});
SH = SH.draw;
SH = SH.legend({'OptRes','PropRes50','PropRes75','PropRes100'});
%xlabel('Arrival rate of slice eMBB (\cdot Default value) [Reqs/sec]');
xlabel('Comp/Comm complexity of slice eMBB (\cdot Default value) [instr/req, bytes/req]');
%xlabel('Routing probability of slice EMBB');
ylabel('Reliability (EMBB)');
%legend('OptRes','PropRes50','PropRes75','PropRes100','Location','best');
%set(gca,'XTickLabel',{'1','10','20','30','40','50','60','70','80','90','100','200','300','400','500'});
figure;
%h3 = bar(scale_factor,prob_success_urllc);
SH = shadowHist(prob_success_urllc,'ShadowType',{'/','.','x','|'});
SH = SH.draw;
SH = SH.legend({'OptRes','PropRes50','PropRes75','PropRes100'});
%xlabel('Arrival rate of slice URLLC (\cdot Default value) [Reqs/sec]');
xlabel('Comp/Comm complexity of slice URLLC (\cdot Default value) [instr/req, bytes/req]');
%xlabel('Routing probability of slice URLLC');
ylabel('Reliability (URLLC)');
%legend('OptRes','PropRes50','PropRes75','PropRes100','Location','best');
%set(gca,'XTickLabel',{'1','10','20','30','40','50','60','70','80','90','100','200','300','400','500'});
xconf = [scale_factor scale_factor(end:-1:1)] ;
yconf = [errpos errneg(end:-1:1)];
yneg = power_opt_final'-errneg;
ypos = errpos - power_opt_final';
figure
%bar(y2,'stacked');
%plot(scale_factor,power_opt,'-k')

hold on;
%p = fill(xconf,yconf,'red');
%p.FaceColor = [1 0.8 0.8]; 
%p.EdgeColor = 'none';
errorbar(scale_factor,power_opt_final,yneg,ypos,"LineStyle","none");
plot(scale_factor,power_opt_final,'r*');
plot(scale_factor,power_prop_50_final,'--k',scale_factor,power_prop_75_final,':k',scale_factor,power_prop_100_final,scale_factor,power_sta_final);


%xlabel('Arrival rate of slices (\cdot Default value) [Reqs/sec]');
xlabel('Comp/Comm complexity of slices (\cdot Default value) [instr/req, bytes/req]');
%xlabel('Routing probability of slices');
ylabel('Total power consumption [Watts]');
legend('95% Confidence Interval','OptRes','PropRes50','PropRes75','PropRes100','Minres','Location','best');
%hold off;
% e.Marker = '*';
% e.MarkerSize = 10;
% e.Color = 'red';

figure;
h3 = bar(power_slice_final,'stacked');
xlabel('Comp/Comm complexity of slices (\cdot Default value) [instr/req, bytes/req]');
ylabel('Power consumption for each slice');
legend('MMTc','eMBB','URLLC','Location','best');