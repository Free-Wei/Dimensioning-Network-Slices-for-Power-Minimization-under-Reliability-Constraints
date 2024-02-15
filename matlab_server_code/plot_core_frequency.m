%dl  = [1, 0.05, 0.005];
%delta = 0.95;
test = [7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140];
net = "Twaren.graphml";
%test = [7,14];
%test = 42;
[G, EdgeTable, NodeTable] = read_graph(net);
%Cpu = ((b-a).*rand(size(NodeTable.Name,1),1)+a)*1e6;  1e12 instructions/sec
Cpu = ones(size(NodeTable.Name,1),1)*1.032e4; 
% 10GHz
Bandwidth = ones(size(EdgeTable.EndNodes,1),1)*1e4;  
%Bandwidth = ((b-a).*rand(size(EdgeTable.EndNodes,1),1)+a)*1e4;
NodeTable = addvars(NodeTable,Cpu);
EdgeTable = addvars(EdgeTable,Bandwidth);
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
    for net_number = 1
    for seed_chose = 10
    power_opt_stochastic = [];
    prob_opt_stochastic = [];
    seed = test(seed_chose);
    for iteration = 1:5
    filename = ['final_duplicated_prob_1s_50ms_5ms_slices_final_',convertStringsToChars(net(net_number)),'_',num2str(seed),'_complexity_delay1e4_iter200_contolerance1-2_inter_10_test_greedy_',num2str(iteration),'.mat'];
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
    filename = ['final_duplicated_prob_1s_50ms_5ms_slices_final_',convertStringsToChars(net(net_number)),'_',num2str(seed),'_complexity_delay1e4_iter200_contolerance1-2_inter_10_test_greedy_',num2str(pos),'.mat'];
    load(filename);
    %%% propotional strategy:
    %[f_prop_50_new{seed}(iter,:),f_prop_75_new{seed}(iter,:),f_prop_100_new{seed}(iter,:),prob_prop_50_new{seed}(iter,:),prob_prop_75_new{seed}(iter,:),prob_prop_100_new{seed}(iter,:)] = load_propstra_com_final(seed,iter);
    %f_prop_50{7}(1,:)
    %prob_prop_50{seed}(iter,:) = prob_prop_50_new{seed}(iter,:);
    %prob_prop_75{seed}(iter,:) = prob_prop_75_new{seed}(iter,:);
    %prob_prop_100{seed}(iter,:) = prob_prop_100_new{seed}(iter,:);
    end
    end
    plot_matrix = zeros(size(NodeTable.Name,1),num_slices);
    for s = 1:num_slices
        for i = 1:size(NodeTable.Name,1)
            plot_matrix(i,s) = x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/3e3;
        end
    end

figure;
x_axis = 1:20;
s = scatter(x_axis,plot_matrix(:,1),'o','CData',[0,0,0]./255,'MarkerFaceAlpha',.5);
hold on;
s1 = scatter(x_axis,plot_matrix(:,2),"+",'CData',[0,0,255]./255,'MarkerFaceAlpha',.5);
hold on;
s2 = scatter(x_axis,plot_matrix(:,3),"square",'CData',[255,0,255]./255,'MarkerFaceAlpha',.5);
ylim([1.5,3.35])
s.SizeData = 150;
s1.SizeData = 150;
s2.SizeData = 150;
xlabel('Physical nodes of physical network "Twaren.graphml"');
ylabel('CPU core frequncy (GHz)');
legend('mMTC','eMBB','URLLC','Location','best');