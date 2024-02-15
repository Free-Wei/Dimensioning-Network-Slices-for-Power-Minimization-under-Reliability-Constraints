%dl  = [1, 0.05, 0.005];
%delta = 0.95;
test = [7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140];
%test = [7,14];
%test = 42;
prob_success_mmtc = zeros(10,4);
prob_success_embb = zeros(10,4);
prob_success_urllc = zeros(10,4);
power_opt = zeros(10,length(test));
power_sta = zeros(10,length(test));
power_prop_50 = zeros(10,length(test));
power_prop_75 = zeros(10,length(test));
power_prop_100= zeros(10,length(test));
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

    for iter = 1:10
    for seed_chose = 1:length(test)
    power_opt_stochastic = [];
    prob_opt_stochastic = [];
    seed = test(seed_chose);
    for iteration = 1:5
    filename = ['final_arrivalrate_1s_50ms_5ms_slices_final_Twaren.graphml_',num2str(seed),'_complexity_delay1e4_iter200_contolerance1-2_inter_',num2str(iter),'_test_greedy_',num2str(iteration),'.mat'];
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
    filename = ['final_arrivalrate_1s_50ms_5ms_slices_final_Twaren.graphml_',num2str(seed),'_complexity_delay1e4_iter200_contolerance1-2_inter_',num2str(iter),'_test_greedy_',num2str(pos),'.mat'];
    load(filename);
    %%% propotional strategy:
    %[f_prop_50_new{seed}(iter,:),f_prop_75_new{seed}(iter,:),f_prop_100_new{seed}(iter,:),prob_prop_50_new{seed}(iter,:),prob_prop_75_new{seed}(iter,:),prob_prop_100_new{seed}(iter,:)] = load_propstra_com_final(seed,iter);
    %f_prop_50{7}(1,:)
    %prob_prop_50{seed}(iter,:) = prob_prop_50_new{seed}(iter,:);
    %prob_prop_75{seed}(iter,:) = prob_prop_75_new{seed}(iter,:);
    %prob_prop_100{seed}(iter,:) = prob_prop_100_new{seed}(iter,:);
    t{iter}(seed_chose,:) = t_app;
    prob_opt_matrix{seed} = (prob_opt{seed}>=0) & (prob_opt{seed}<=1.0001);
    prob_success_mmtc(iter,1) = prob_success_mmtc(iter,1) + prob_opt{seed}(iter,1)*prob_opt_matrix{seed}(iter,1);
    prob_success_embb(iter,1) = prob_success_embb(iter,1) + prob_opt{seed}(iter,2)*prob_opt_matrix{seed}(iter,2);
    prob_success_urllc(iter,1) = prob_success_urllc(iter,1) + prob_opt{seed}(iter,3)*prob_opt_matrix{seed}(iter,3);
    
    prob_50{seed} = (prob_prop_50_new{seed}>=0) & (prob_prop_50_new{seed}<=1.0001);
    if prob_50{seed}(iter,1) == 1
    prob_success_mmtc(iter,2) = prob_success_mmtc(iter,2) + prob_prop_50_new{seed}(iter,1)*prob_50{seed}(iter,1);
    end
    if prob_50{seed}(iter,2) == 1
    prob_success_embb(iter,2) = prob_success_embb(iter,2) + prob_prop_50_new{seed}(iter,2)*prob_50{seed}(iter,2);
    end
    if prob_50{seed}(iter,3) == 1
    prob_success_urllc(iter,2) = prob_success_urllc(iter,2) + prob_prop_50_new{seed}(iter,3)*prob_50{seed}(iter,3);
    end
    prob_75{seed} = (prob_prop_75_new{seed}>=0) & (prob_prop_75_new{seed}<=1.0001);
    if prob_75{seed}(iter,1) == 1
    prob_success_mmtc(iter,3) = prob_success_mmtc(iter,3) + prob_prop_75_new{seed}(iter,1)*prob_75{seed}(iter,1);
    end
    if prob_75{seed}(iter,2) == 1
    prob_success_embb(iter,3) = prob_success_embb(iter,3) + prob_prop_75_new{seed}(iter,2)*prob_75{seed}(iter,2);
    end
    if prob_75{seed}(iter,3) == 1
    prob_success_urllc(iter,3) = prob_success_urllc(iter,3) + prob_prop_75_new{seed}(iter,3)*prob_75{seed}(iter,3);
    end
    
    prob_100{seed} = (prob_prop_100_new{seed}>=0) & (prob_prop_100_new{seed}<=1.0001);
    if prob_100{seed}(iter,1) == 1
    prob_success_mmtc(iter,4) = prob_success_mmtc(iter,4) + prob_prop_100_new{seed}(iter,1)*prob_100{seed}(iter,1);
    end
    if prob_100{seed}(iter,2) == 1
    prob_success_embb(iter,4) = prob_success_embb(iter,4) + prob_prop_100_new{seed}(iter,2)*prob_100{seed}(iter,2);
    end
    if prob_100{seed}(iter,3) == 1
    prob_success_urllc(iter,4) = prob_success_urllc(iter,4) + prob_prop_100_new{seed}(iter,3)*prob_100{seed}(iter,3);
    end
    power_opt(iter,seed_chose) = sum(f_opt{seed}(iter,:));
    power_sta(iter,seed_chose) = sum(f_sta{seed}(iter,:));
    power_prop_50(iter,seed_chose) = sum(f_prop_50{seed}(iter,:));  
    power_prop_75(iter,seed_chose) = sum(f_prop_75{seed}(iter,:));  
    power_prop_100(iter,seed_chose) = sum(f_prop_100{seed}(iter,:));  
        for n = 1:num_slices
            power_slice_final(iter,n) = power_slice_final(iter,n) + f_opt{seed}(iter, n);
        end
    end
    prob_success_mmtc(iter,:) = prob_success_mmtc(iter,:)/length(test);
    prob_success_embb(iter,:) = prob_success_embb(iter,:)/length(test);
    prob_success_urllc(iter,:) = prob_success_urllc(iter,:)/length(test);
    power_opt_final(iter) = mean(power_opt(iter,:));
    power_sta_final(iter) = mean(power_sta(iter,:));
    power_prop_50_final(iter) = mean(power_prop_50(iter,:));
    power_prop_75_final(iter) = mean(power_prop_75(iter,:));
    power_prop_100_final(iter) = mean(power_prop_100(iter,:));
    [muhat,sigmahat,muci,sigmaci] = normfit(power_opt(iter,:));
    errneg(iter) = muci(1);
    errpos(iter) = muci(2);

    end
    power_slice_final =   power_slice_final /length(test);
fprintf( 'drop percentage: %f \n', count_drop/count_all);
y = [prob_success_mmtc(1,:);prob_success_embb(1,:);prob_success_urllc(1,:)];
%y2 = [power_sta,power_opt,power_prop_50,power_prop_100,]';

%legend('OptRes','PropRes50','PropRes75','PropRes100','Location','best');
%set(gca,'XTickLabel',{'1','10','20','30','40','50','60','70','80','90','100','200','300','400','500'});
xconf = [scale_factor scale_factor(end:-1:1)] ;
yconf = [errpos errneg(end:-1:1)];
figure;
%bar(y2,'stacked');
%plot(scale_factor,power_opt,'-k')

hold on;
p = fill(xconf,yconf,'red');
p.FaceColor = [1 0.8 0.8]; 
p.EdgeColor = 'none';

plot(scale_factor,power_opt_final,'r*');
plot(scale_factor,power_prop_50_final,'--k',scale_factor,power_prop_75_final,':k',scale_factor,power_prop_100_final,scale_factor,power_sta_final);

% e = errorbar(scale_factor,power_opt_final,err);
xlabel('Ingress arrival rate of slices (\cdot Default value) [Reqs/sec]');
%xlabel('Comp/Comm complexity of slices (\cdot Default value) [instr/req, bytes/req]');
%xlabel('Routing probability of slices');
ylabel('Overall power consumption [Watts]');
legend('95% Confidence Interval','OptRes','PropRes50','PropRes75','PropRes100','Minres','Location','best');
%hold off;
% e.Marker = '*';
% e.MarkerSize = 10;
% e.Color = 'red';
matrix2draw2 = prob_success_mmtc;
threshold_value = 0.9;
matrix2draw2(matrix2draw2<threshold_value) = threshold_value;
x_range = 1:10;
%y_range = repmat(threshold_value,10,1)';
figure;
%plot(x_range,y_range,'LineWidth',2);
%hold on;
%shadowFill(x_range,y_range,pi/4,80)
plot(x_range-0.2, matrix2draw2(:,1)-0.05,'o')
hold on;
plot(x_range-0.1,matrix2draw2(:,2)-0.05,'+')
hold on;
plot(x_range,matrix2draw2(:,3)-0.05,'x')
hold on;
plot(x_range+0.1,matrix2draw2(:,4)-0.05,'pentagram')
hold on;
yline(threshold_value,'-','Realibility')
ylim([threshold_value-0.1 1.0])
%xticks([1,2,3,4,5,6,7,8,9,10,11]);
%set(gca,'XTickLabel',{'0.3','0.35','0.4','0.45','0.5','0.55','0.6',"0.65","0.7","0.75"})
%xticklabels({'0.3','0.35','0.4','0.45','0.5','0.55','0.6',"0.65","0.7","0.75"})
yt = get(gca, 'YTick');
threshold_axis = find(yt==threshold_value);
set(gca, 'YTick',yt(threshold_axis:end), 'YTickLabel',sprintfc('%.2f',yt(threshold_axis:end)))
xlabel('Ingress arrival rate of slices (\cdot 5000 + 5000) [Reqs/sec]');
%xlabel('Routing probability of slices');
ylabel('Success or failure');
legend('OptRes','PropRes50','PropRes75','PropRes100');


matrix2draw1 = prob_success_embb;
threshold_value = 0.95;
matrix2draw1(matrix2draw1<threshold_value) = threshold_value;
x_range = 1:10;
y_range = repmat(threshold_value,10,1)';
figure;
%plot(x_range,y_range,'LineWidth',2);
%hold on;
%shadowFill(x_range,y_range,pi/4,80);
plot(x_range-0.2, matrix2draw1(:,1)-0.025,'o')
hold on;
plot(x_range-0.1,matrix2draw1(:,2)-0.025,'+')
hold on;
plot(x_range,matrix2draw1(:,3)-0.025,'x')
hold on;
plot(x_range+0.1,matrix2draw1(:,4)-0.025,'pentagram')
hold on;
yline(threshold_value,'-','Realibility')
ylim([threshold_value-0.05 1.0])
yt = get(gca, 'YTick');
threshold_axis = find(yt==threshold_value);
set(gca, 'YTick',yt(threshold_axis:end), 'YTickLabel',sprintfc('%.2f',yt(threshold_axis:end)))
xlabel('Ingress arrival rate of slices (\cdot 500) [Reqs/sec]');
%xlabel('Routing probability of slices');
ylabel('Success or failure');
legend('OptRes','PropRes50','PropRes75','PropRes100');

matrix2draw = prob_success_urllc;
threshold_value = 0.99;
matrix2draw(matrix2draw<threshold_value) = threshold_value;
x_range = 1:10;
y_range = repmat(threshold_value,10,1)';
figure;
%plot(x_range,y_range,'LineWidth',2);
%hold on;
%shadowFill(x_range,y_range,pi/4,80);
plot(x_range-0.2, matrix2draw(:,1)-0.005,'o')
hold on;
plot(x_range-0.1,matrix2draw(:,2)-0.005,'+')
hold on;
plot(x_range,matrix2draw(:,3)-0.005,'x')
hold on;
plot(x_range+0.1,matrix2draw(:,4)-0.005,'pentagram')
hold on;
yline(threshold_value,'-','Realibility')
ylim([threshold_value-0.01 1.0])
yt = get(gca, 'YTick');
threshold_axis = find(yt==threshold_value);
set(gca, 'YTick',yt(threshold_axis:end), 'YTickLabel',sprintfc('%.2f',yt(threshold_axis:end)))
xlabel('Ingress arrival rate of slices (\cdot 1000 + 3000) [Reqs/sec]');
%xlabel('Routing probability of slices');
ylabel('Success or failure');
legend('OptRes','PropRes50','PropRes75','PropRes100');


figure;

%SH = shadowHist(power_slice_final,'stacked','ShadowType',{'/','.','x'});
%SH = SH.draw;
%SH = SH.legend({'OptRes','PropRes50','PropRes75','PropRes100'});
%xlabel('Arrival rate of slice mMTC (\cdot Default value) [Reqs/sec]');
%xlabel('Comp/Comm complexity of slice mMTC (\cdot Default value) [instr/req, bytes/req]');
%xlabel('Routing probability of slice mMTC');
%ylabel('Reliability (mMTC)');
%legend('OptRes','PropRes50','PropRes75','PropRes100','Location','best');
h3 = bar(power_slice_final,'stacked');
xlabel('Ingress requests rates of slice');
ylabel('Power consumption for each slice');
legend('MMTc','eMBB','URLLC','Location','best');
