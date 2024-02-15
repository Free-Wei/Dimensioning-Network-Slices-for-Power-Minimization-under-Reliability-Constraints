%read graph and define capacity for links and nodes
[G, EdgeTable, NodeTable] = read_graph("Twaren.graphml");
a = 1;
b = 5;
%Cpu = ((b-a).*rand(size(NodeTable.Name,1),1)+a)*1e6;  1e12 instructions/sec
Cpu = ones(size(NodeTable.Name,1),1)*1.032e4; 
% 10GHz
Bandwidth = ones(size(EdgeTable.EndNodes,1),1)*1e4;  
%Bandwidth = ((b-a).*rand(size(EdgeTable.EndNodes,1),1)+a)*1e4;
NodeTable = addvars(NodeTable,Cpu);
EdgeTable = addvars(EdgeTable,Bandwidth);
f_sta ={};
f_prop_50= {};
f_prop_75 = {};
f_prop_100 = {};
f_opt = {};
T_prop_50 = {};
T_prop_75 ={};
T_prop_100 = {};
T_opt = {};
prob_prop_50 = {};
prob_prop_75 ={};
prob_prop_100 ={};
prob_opt = {};
rep ={};
base = 0;
%test = [7,14,21,28,35,43,58,63,72,78,84,91,98,105,112,119,126,133,140];
test = 126;
%test = [7,14,21,28,35,42,49,56,63,70];
%test = 56;
%test = [98,105,112,119,126,133,140];
%test = [112,119,126,13 143,140];
%test = 140;
%test different seeds
%for seed = 7
%remake = randi([5 15],3,1);
%remake
for seed_chose = 1:length(test)
seed = test(seed_chose);
seed = seed_verify(seed);
rng(seed);

%flag_overload = 1;
%base = 0;
%while(flag_overload == 1)
%     base = base + 1;
%     rng(seed+base);
% fprintf( 'iteration times: %d \n', base);
num_slices = 3;
min_num_components = 2;
max_num_components = 5;
max_number_of_replicas = 2;
max_number_of_ingress = 2;
num_components_slice = [randi([min_num_components,max_num_components],1,num_slices-1),randi([2,3],1,1)];
dl = [1,0.05,0.005]*1e4;

%for num_iteration = 1:10
replicas_destinations = {};
probability = {};
rng(seed);
for s = 1:num_slices
    replicas_destinations{s} = {};
    probability_s = [];
    len_components= [];
    listForSample = (1:size(NodeTable.Name,1));
    %define the ingress nodes
    w_dup = 0.5;
    %rng(seed);
    replicas_destinations{s}{1} = randsample(listForSample, randsample(max_number_of_ingress,1,true,[0 1]));
    if size(replicas_destinations{s}{1},2) > 1
        probability_s(1) = 1;
    end
    for i = 2:num_components_slice(s)
        p = 0;
        replicas_destinations{s}{i} = randsample(listForSample, randsample(max_number_of_replicas,1,true,[1-w_dup w_dup] ));
        %replicas_destinations{s}{i} = randsample(listForSample, 1);
        len_components(s) = length(replicas_destinations{s}{i});
        for l = 1:size(replicas_destinations{s}{i},2)
            val = find(listForSample == replicas_destinations{s}{i}(l));
            listForSample(val) = [];
        end
        if size(replicas_destinations{s}{i},2) > 1
            p = 1;
        end
        probability_s(i) = p;
    end
    probability{s} = probability_s;
end
%replicas_destinations{2}
%end

%find the shortestpath between 2 nodes;
paths = {};
for i1 = 1: size(NodeTable.Name,1)
   paths{i1} = {};
   for i2 = 1: size(NodeTable.Name,1)
       if i1 ~= i2
           paths{i1}{i2} = shortestpath(G, i1, i2);
       else
           paths{i1}{i2} = [];
       end
   end
end

%traffic for nodes
%arrival_rate = (10:(30-10)/(num_slices-1):30);
for num_iteration = 1:11

coefficient_bw1 = 14.55/550;
coefficient_bw2 = 4.5;
coefficient_bw3 = 19.055;
beta_max = [1e-3,1.2e-2,2.4e-3];
alpha_max = [6e-1, 7.2, 1.44];
beta_list = [1.6e-4, 5.12e-4, 2.56e-4 ];
alpha_list = [9.6e-2, 3.072e-1, 1.536e-1];
%alpha_list = [9.5e-2, 5.5296, 1.0752];
diffrent_beta = (beta_max-beta_list)/10;
diffrent_alpha = (alpha_max-alpha_list)/100;
%coefficient_alpha = rand(size(NodeTable.Name,1),1)*3e2;
% 3e8 instructions/ request
%coefficient_alpha = ones(size(NodeTable.Name,1),1)*3e2;
%flag = [500 250 150];
for s = 1:num_slices
    coefficient_alpha(:,s) = ones(size(NodeTable.Name,1),1)*(alpha_list(s)+diffrent_alpha(s)*(20+5*5));
    %size(coefficient_beta(:,s))
    coefficient_alpha(1,s)
    %size(ones(size(EdgeTable.EndNodes,1),1)*beta_list(s))
    coefficient_beta(:,s) = ones(size(EdgeTable.EndNodes,1),1)*(beta_list(s)+diffrent_beta(s)*5);
    coefficient_beta(1,s)
end
%coefficient_beta = rand(size(EdgeTable.EndNodes,1),1)*8.5e-2;
% 85Kb/request
%coefficient_beta = ones(size(EdgeTable.EndNodes,1),1)*8.5e-2;
for greedy_num = 1:5
%num_iteration = 5;
prob_fail = [0.1,0.05,0.01];
%rate = [100];
%for r = 1:length(rate)
%arrival_rate = [];
%for s = 1:num_slices
%for multi_value = 1:9
%arrival_rate = [100*multi_value 20*multi_value 5*multi_value];
%arrival_rate = [100 20 5]*num_iteration;
%arrival_rate = [10000 1000 5000] + [500 40 200]*5;
arrival_rate = [5000 100 3000]+ [1000 80 200]*5; 
%p1 = [0.5,0.5];
%p_route = 0.5;
p_route = 0.1*(num_iteration-1);
prob_routing = [p_route 1-p_route];
traffic_path = {};
traffic_arrival = {};
temp_traffic = {};
traffic_path_rep ={};
for i = 1:num_slices
    traffic_path_rep{i} = zeros(max_number_of_replicas,size(EdgeTable.EndNodes,1)+size(NodeTable.Name,1));
    total_arrival(i) = size(replicas_destinations{i}{1},2)*arrival_rate(i);
    temp_traffic{i} = zeros(2,num_components_slice(i));
end 
routing_paths ={};
prob_paths = {};
for s = 1:num_slices  
    routing_paths{s} = {};
    prob_paths{s} = {};
    %flag = cell(1,num_components_slice(s)-1);
    node_decide{s} = zeros(size(NodeTable.Name,1),size(NodeTable.Name,1));
    for i = 1:num_components_slice(s)-1
        %flag{i} = zeros(1,size(EdgeTable.EndNodes,1));
        start_point = replicas_destinations{s}{i};
        finish_point = replicas_destinations{s}{i+1};
        for start = 1:size(start_point,2)
            path_length = ones(1,size(NodeTable.Name,1))*1e5;
            for finish = 1:size(finish_point,2)
                path_length(finish) = length(paths{start_point(start)}{finish_point(finish)});
            end
            if length(unique(path_length)) == 2 && size(finish_point,2) == 2
                for finish = 1:size(finish_point,2)
                    %node_decide{s}(start_point(start),finish_point(finish)) = 0.5;
                    node_decide{s}(start_point(start),finish_point(finish)) = prob_routing(finish);
                    if finish == 1
                    routing_paths{s}{i}{start} = paths{start_point(start)}{finish_point(finish)};
                    prob_paths{s}{i}{start} = prob_routing(finish);
                    else
                    %node_decide{s}(start_point(start),finish_point(finish)) = prob_routing(finish);
                    routing_paths{s}{i}{start} = [routing_paths{s}{i}{start}; paths{start_point(start)}{finish_point(finish)}];
                    prob_paths{s}{i}{start} = [prob_paths{s}{i}{start};prob_routing(finish)];
                    end
                end
            else
                [q,w] = min(path_length);
                 node_decide{s}(start_point(start),finish_point(w)) = 1;
                 if isempty(paths{start_point(start)}{finish_point(w)}) == 1
                        routing_paths{s}{i}{start} = start_point(start);
                        prob_paths{s}{i}{start} = 1;
                 else
                    routing_paths{s}{i}{start}=  paths{start_point(start)}{finish_point(w)};
                    prob_paths{s}{i}{start} = 1;
                 end
            end
        end
    end
end
%%
for s = 1:num_slices
    for i = 1:num_components_slice(s)
        for l = 1:size(replicas_destinations{s}{i},2)
            if i == 1
               traffic_arrival{s}(l,replicas_destinations{s}{i}(l)) = arrival_rate(s);    
            elseif probability{s}(i-1) == 1 && probability{s}(i) == 1
                destination_decide = sum(node_decide{s},1);
                for j = 1:size(replicas_destinations{s}{i-1},2)
                    if node_decide{s}(replicas_destinations{s}{i-1}(j),replicas_destinations{s}{i}(l)) ~= 0 
                        %if  (ismember(replicas_destinations{s}{i}(l),replicas_destinations{s}{1})  == 1 && (i == 2))
                        
                    %temp_val = sum(temp_traffic{s},1);
                        %traffic_path{s}(l,replicas_destinations{s}{i}(l)) = temp_val(i-1)*p1(l);
                            %if (replicas_destinations{s}{i-1}(j) ~= replicas_destinations{s}{i}(l))
                            %    traffic_path{s}(l,replicas_destinations{s}{i}(l)) =  traffic_path{s}(l,replicas_destinations{s}{i}(l)) + temp_traffic{s}(j,i-1)*node_decide{s}(replicas_destinations{s}{i-1}(j),replicas_destinations{s}{i}(l));
                            %end
                        %elseif (ismember(replicas_destinations{s}{i}(l),replicas_destinations{s}{1})  == 1 && (i > 2))
                        %    traffic_path{s}(l,replicas_destinations{s}{i}(l)) = temp_traffic{s}(j,i-1)*node_decide{s}(replicas_destinations{s}{i-1}(j),replicas_destinations{s}{i}(l));
                        %else
                            traffic_path_rep{s}(l,replicas_destinations{s}{i}(l)) = traffic_path_rep{s}(l,replicas_destinations{s}{i}(l)) + temp_traffic{s}(j,i-1)*node_decide{s}(replicas_destinations{s}{i-1}(j),replicas_destinations{s}{i}(l));
                        %end
                    end
                end
            elseif probability{s}(i-1) == 0 && probability{s}(i) == 1
                    if node_decide{s}(replicas_destinations{s}{i-1}(1),replicas_destinations{s}{i}(l)) ~= 0
                        traffic_path_rep{s}(l,replicas_destinations{s}{i}(l)) = temp_traffic{s}(1,i-1)*node_decide{s}(replicas_destinations{s}{i-1}(1),replicas_destinations{s}{i}(l));
                    end
            elseif probability{s}(i-1) == 1 && probability{s}(i) == 0
%                     if traffic_path_rep{s}(l+1,replicas_destinations{s}{i-1}(l+1)) ~= 0
%                         traffic_path_rep{s}(l+1,replicas_destinations{s}{i-1}(l+1)) = 0;
%                     end
                    temp_val = sum(temp_traffic{s},1);
                    traffic_path_rep{s}(l,replicas_destinations{s}{i}(l)) = temp_val(i-1);
            else
               traffic_path_rep{s}(l,replicas_destinations{s}{i}(l)) = temp_traffic{s}(l,i-1);
            end
        end
        for n = 1:size(replicas_destinations{s}{i},2)
            if i == 1
               temp_traffic{s}(n,i) = traffic_arrival{s}(n,replicas_destinations{s}{i}(n));
            else
               temp_traffic{s}(n,i) = traffic_path_rep{s}(n,replicas_destinations{s}{i}(n));
            end
        end
    end
    traffic_path{s} = sum(traffic_path_rep{s},1);
end
%%
%traffic for links:
for s = 1:num_slices  
    [G, EdgeTable, NodeTable] = read_graph("Twaren.graphml");
    NodeTable = addvars(NodeTable,Cpu);
    EdgeTable = addvars(EdgeTable,Bandwidth);
    %flag = cell(1,num_components_slice(s)-1);
    for i = 1:num_components_slice(s)-1
        %flag{i} = zeros(1,size(EdgeTable.EndNodes,1));
        start_point = replicas_destinations{s}{i};
        finish_point = replicas_destinations{s}{i+1};
        for start = 1:size(start_point,2)
            for finish = 1:size(finish_point,2)
                if size(finish_point,2) == 1
                    %paths{i1}{i2} = shortestpath(G, start_point(start), finish_point(finish)};);
                    %path_between_them = paths{start_point(start)}{finish_point(finish)};
                    path_between_them =  shortestpath(G, start_point(start), finish_point(finish));
                    if i == 1
                        traffic = arrival_rate(s);
                    else
                        traffic = traffic_path{s}(start_point(start));
                    end
                    for i_path = 1:size(path_between_them,2)-1
                        val = find(EdgeTable.EndNodes(:,1) == path_between_them(i_path) & EdgeTable.EndNodes(:,2) == path_between_them(i_path+1));
                        %val_G = find(G.Edges.EndNodes(:,1) == num2str(path_between_them(i_path)) && G.Edges.EndNodes(:,2) == path_between_them(i_path+1));
                        %G = rmedge(G, path_between_them(i_path),path_between_them(i_path+1));
                        %flag{i}(val) = flag{i}(val) + traffic;
                        %if traffic_path{s}(val+size(NodeTable.Name,1)) == 0 
                        traffic_path{s}(val+size(NodeTable.Name,1)) = traffic_path{s}(val+size(NodeTable.Name,1)) + traffic;
                        %traffic_path{s}(val+size(NodeTable.Name,1)) = traffic;
                        %flag_overload = 0;
                        %else
                        %   flag_overload = 1;
                        %    break;
                        %end
                    end
                elseif size(start_point,2) == 1 && size(finish_point,2) == 2
                    %path_between_them = paths{start_point(start)}{finish_point(finish)};
                    path_between_them =  shortestpath(G, start_point(start), finish_point(finish));
                    traffic = traffic_path{s}(finish_point(finish));
                    for i_path = 1:size(path_between_them,2)-1
                        val = find(EdgeTable.EndNodes(:,1) == path_between_them(i_path) & EdgeTable.EndNodes(:,2) == path_between_them(i_path+1));
                        %flag{i}(val) = flag{i}(val) + traffic;
                        %if traffic_path{s}(val+size(NodeTable.Name,1)) == 0
                        %val_G = find(G.Edges.EndNodes(:,1) == path_between_them(i_path) & G.Edges.EndNodes(:,2) == path_between_them(i_path+1));
                       
                        traffic_path{s}(val+size(NodeTable.Name,1)) = traffic_path{s}(val+size(NodeTable.Name,1)) + traffic;
                        %if traffic_path{s}(val+size(NodeTable.Name,1)) > 0
                        %   G = rmedge(G, path_between_them(i_path),path_between_them(i_path+1));
                        %end
                        %traffic_path{s}(val+size(NodeTable.Name,1)) = traffic;
                        %flag_overload = 0;
                        %else
                        %    flag_overload = 1;
                        %    break;
                        %end
                    end
                elseif node_decide{s}(start_point(start),finish_point(finish)) ~= 0
                    %path_between_them = paths{start_point(start)}{finish_point(finish)};
                    path_between_them = shortestpath(G, start_point(start), finish_point(finish));
                        if i == 1
                            traffic = arrival_rate(s)*node_decide{s}(start_point(start),finish_point(finish));
                        else
                            traffic = traffic_path{s}(start_point(start))*node_decide{s}(start_point(start),finish_point(finish));
                        end
                    for i_path = 1:size(path_between_them,2)-1
                        val = find(EdgeTable.EndNodes(:,1) == path_between_them(i_path) & EdgeTable.EndNodes(:,2) == path_between_them(i_path+1));
                        %flag{i}(val) = flag{i}(val) + traffic
                        traffic_path{s}(val+size(NodeTable.Name,1)) = traffic_path{s}(val+size(NodeTable.Name,1)) + traffic;
                        %if traffic_path{s}(val+size(NodeTable.Name,1)) == 0
                        %%traffic_path{s}(val+size(NodeTable.Name,1)) = traffic;
                        %%%%flag_overload = 0;
                        %else
                        %    flag_overload = 1;
                        %    break;
                        %end
                    end                    
                end
            end
        end
        for start = 1:size(start_point,2)
            for finish = 1:size(finish_point,2)
                path_between_them =  shortestpath(G, start_point(start), finish_point(finish));
                for i_path = 1:size(path_between_them,2)-1
                    val = find(EdgeTable.EndNodes(:,1) == path_between_them(i_path) & EdgeTable.EndNodes(:,2) == path_between_them(i_path+1));
                    if traffic_path{s}(val+size(NodeTable.Name,1)) > 0
                        G = rmedge(G, path_between_them(i_path),path_between_them(i_path+1));
                    end
                end
            end
        end
    end
%     %for j = 1:num_components_slice(s)-1
%     %    for v = 1:size(EdgeTable.EndNodes,1)
%             %traffic_path{s}(v+size(NodeTable.Name,1)) = max(traffic_path{s}(v+size(NodeTable.Name,1)),flag{j}(v));
%             if traffic_path{s}(v+size(NodeTable.Name,1)) > total_arrival(s)
%                 flag_overload = 1;
%                 break;
%             end
%         end
%     end
end
%end




% capacity for links and nodes
Coeff = [zeros(size(NodeTable.Name,1),1);EdgeTable.Bandwidth]';
for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)
        if traffic_path{s}(i) == 0 || (ismember(i,replicas_destinations{s}{1})&&(~ismember(i,[replicas_destinations{s}{2:end}])))
            x_range_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
            %x_sta(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            x_range_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = ((coefficient_alpha(i,s)*traffic_path{s}(i)));
            %x_sta(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (coefficient_alpha(i)*traffic_path{s}(i))+1;
        end
    end
    for k = 1:size(EdgeTable.EndNodes,1)
        if traffic_path{s}(k+size(NodeTable.Name,1)) == 0
            x_range_min(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
            %x_sta(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            x_range_min(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (coefficient_beta(i,s)*traffic_path{s}(k+size(NodeTable.Name,1)))+1e-6;
            %x_sta(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (coefficient_beta(i)*traffic_path{s}(k+size(NodeTable.Name,1)))+1;
        end
    end
end
x_range_max = [];
for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)
        if  traffic_path{s}(i) == 0 || (ismember(i,replicas_destinations{s}{1})&&(~ismember(i,[replicas_destinations{s}{2:end}])))
            x_range_max(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            x_range_max(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = NodeTable.Cpu(i);
        end
    end
    for k = 1:size(EdgeTable.EndNodes,1)
        if traffic_path{s}(k+size(NodeTable.Name,1)) == 0
            x_range_max(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            x_range_max(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = EdgeTable.Bandwidth(k);
        end
    end
end
%x_range_max = [Coeff Coeff Coeff Coeff Coeff];

% constarints like x1 + x61 + x121 <= RV
B = diag([zeros(1,size(NodeTable.Name,1)),ones(1,size(EdgeTable.EndNodes,1))]);
A = [];
for s = 1:num_slices
    A = [A,B];
end
%A = [B B B B B];
%A = repmat(B,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices,size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1));
b = Coeff';

% initialise x0
for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)
        if  traffic_path{s}(i) == 0 || (ismember(i,replicas_destinations{s}{1})&&(~ismember(i,[replicas_destinations{s}{2:end}])))
           x0(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
           x0(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) =  (coefficient_alpha(i,s)*traffic_path{s}(i))+rand(1);
        end
    end
    for k = 1:size(EdgeTable.EndNodes,1)
       if traffic_path{s}(k+size(NodeTable.Name,1)) == 0
            x0(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
       else
            x0(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (coefficient_beta(i,s)*traffic_path{s}(k+size(NodeTable.Name,1)))+rand(1);
       end
    end
end

%%%%phase 1: r_max = NodeTable.Cpu(i)
%fun2 = @(x)delay_large_hyex(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations);
r_max = [];
for s = 1:num_slices
    r_max = [r_max, [NodeTable.Cpu',zeros(1,size(EdgeTable.EndNodes,1))] ];
end
epsilon = ones(1,num_slices)*1e6;
round = ones(1,num_slices)*1;
step = 0;
threshold = [1e-2,5e-3,1e-3];
unchanged_times = 3;
buffer_check = zeros(num_slices,unchanged_times);
upp_round = ones(1,num_slices)*2;
low_round = ones(1,num_slices);
t_start = tic;
while (sum(abs(epsilon) > threshold ) ~= 0)
rng('shuffle');
eps = rand(1,3);
eps
[t_app,round,step] = approxiamate_t_dvfs(dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations,x_range_max,step,round);
%round
%upp_round
%low_round
%step = step + 1;
fun2 = @(x)delay_hoeffding(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,replicas_destinations,x_range_max,t_app);
opts = optimoptions('fmincon','Display','iter','Algorithm','interior-point','SpecifyConstraintGradient',true);
%opts = optimoptions('fmincon','SpecifyConstraintGradient',true);
opts = optimoptions(opts,'MaxFunctionEvaluations',9e6,'MaxIterations',500, 'StepTolerance',1e-10,'ConstraintTolerance',1e-6); % Recommended

[x,fval,exitflag,output] = fmincon(@(x)ObjFunc(x,traffic_path,coefficient_alpha,NodeTable,EdgeTable,num_slices),x0,A,b,[],[],x_range_min,x_range_max, fun2, opts);
%x = x +rand(1);
[non,res] = delay_hoeffding(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,replicas_destinations,x_range_max,t_app);
%prob = delay_large_hyex(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations);
prob = delay_large_real_hyex(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,routing_paths,prob_paths,replicas_destinations,num_components_slice);
buffer_check(mod(step,unchanged_times)+1,:) = prob;
if(((step>unchanged_times && size(uniquetol(buffer_check,1e-4,'ByRows',true),1) == 1)||step>=50))
    break;
end
%if step >= 100
%    break;
%end
round
step
prob
greedy = 0.07;
step = step + 1;
flag_lower = zeros(1,num_slices);
for n = 1: length(prob)
    epsilon(n) = 1-prob_fail(n) - prob(n);
    if prob(n) < 0 || prob(n) > 1 
        eps(n) = 0.5 *rand(1);
    elseif isnan(prob(n))
        eps(n) = 1e-4;
    end
    if (epsilon(n) < 0 && abs(epsilon(n)) > threshold(n) ) || isnan(epsilon(n))
           if (epsilon(n) < 0 && abs(epsilon(n)) > threshold(n) ) || isnan(epsilon(n))
            if eps(n) <= 0.07 && flag_lower(n) ~= 0
                upp_round(n) = flag_upper(n);
                low_round(n) = 1;
                round(n) = (upp_round(n)-low_round(n))/2+low_round(n);
            elseif  abs(low_round(n)-upp_round(n)) <=1e-4 && flag_lower(n) ~= 0
                upp_round(n) = flag_upper(n);
                round(n) = (upp_round(n)-low_round(n))/2+low_round(n); 
            elseif flag_lower(n) == 0
                upp_round(n) = 2*round(n);
                low_round(n) = round(n);
                round(n) = (upp_round(n)-low_round(n))/2+low_round(n);
                flag_upper(n)= upp_round(n);
            else
                low_round(n) = round(n);
                round(n) = (upp_round(n)-low_round(n))/2+low_round(n); 
            end
    elseif epsilon(n) > 0 && abs(epsilon(n)) > threshold(n)
            flag_lower(n) = flag_lower(n) + 1;
            if eps(n) <= 0.07
                %upp_round(n) = 2*round(n);
                low_round(n) = 1;
                round(n) = (upp_round(n)-low_round(n))/2+low_round(n);
            elseif abs(low_round(n)-upp_round(n)) <=1e-4
                low_round(n) = 1;
                round(n) = (upp_round(n)-low_round(n))/2+low_round(n); 
            %elseif ismember(step, remake(n,:))
            %    low_round(n) = 0;
            %    round(n) = (upp_round(n)-low_round(n))/2+low_round(n); 
            else
            upp_round(n) = round(n);
            round(n) = (upp_round(n)-low_round(n))/2+low_round(n);
            end
           end
    end
end
%upp_round
%low_round
% if exitflag == -2 || exitflag == 0
%     step = step -1;
%     for n = 1:length(prob)
%         round(n) = round(n) + 1e-4;
%     end
% else
%     for n = 1: length(prob)
%         epsilon(n) = 1-prob_fail(n) - prob(n);
%         if epsilon(n) < 0 && abs(epsilon(n)) > threshold(n)
%             if step <= 13
%                 round(n) = round(n)+1/(2^(step));
%             else
%                 round(n) = round(n)+1/(2^(13));
%             end
%         elseif epsilon(n) > 0 && abs(epsilon(n)) > threshold(n)
%             if step <= 13
%                 round(n) = round(n)-1/(2^(step));
%             else
%                 round(n) = round(n)-1/(2^(13));
%             end
%         end
%     end
% end

% if step > 10
%     break;
% end
end
prob
res
fval
% delay_cost = zeros(1,size(res,2));
% for m = 1:size(res,2)
%     delay_cost(m) = res(m)+(dl(m)/(-log(prob_fail(m))));
% end
%%%%%phase 2: change r_max
possible_p = [1.67,1.89,2.09,2.31,2.51,2.70,2.90,3.07,3.26,3.44]*3e3;
for i = 1:num_slices
    for j = 1:size(NodeTable.Name,1)
        if  x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) ~= 0
            index = find(possible_p>x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))));
            x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))  = possible_p(index(1))+rand(1);
        else
            x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        end
    end
end

[non,res] = delay_hoeffding(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,replicas_destinations,x_range_max,t_app);
%prob = delay_large_hyex(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations);
prob = delay_large_real_hyex(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,routing_paths,prob_paths,replicas_destinations,num_components_slice);
prob
res
fval
tEnd = toc(t_start);




%other strategies
x_range_sta = x_range_min;
x_range_prop_50 = 0.5*x_range_max;
x_range_prop_75 = 0.75*x_range_max;
x_range_prop_100 = x_range_max;
count = zeros(num_slices,size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1));
for i = 1:num_slices
    for k = 1:size(EdgeTable.EndNodes,1)
        if  x_range_prop_100(k+size(NodeTable.Name,1)+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) ~= 0
            count(i,k+size(NodeTable.Name,1)) = coefficient_beta(1,i)*traffic_path{i}(k+size(NodeTable.Name,1));
        end
    end
end
for s = 1:num_slices
    for j = 1:size(NodeTable.Name,1)
        x_range_prop_50(j+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0.5*x_range_max(j+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)));
        x_range_prop_75(j+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0.75*x_range_max(j+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)));
        x_range_prop_100(j+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = x_range_max(j+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)));
    end
    for m = size(NodeTable.Name,1):(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))
        if sum(count(:,m)) ~= 0
            x_range_prop_50(m+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = count(s,m) * 0.5*x_range_max(m+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/sum(count(:,m));
            x_range_prop_75(m+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = count(s,m) * 0.75*x_range_max(m+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/sum(count(:,m));
            x_range_prop_100(m+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = count(s,m) * x_range_max(m+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/sum(count(:,m));
        end
    end
end

for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)
        if  traffic_path{s}(i) == 0 || (ismember(i,replicas_destinations{s}{1})&&(~ismember(i,[replicas_destinations{s}{2:end}])))
            x_range_sta(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
            x_range_prop_50(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
            x_range_prop_75(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
            x_range_prop_100(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else 
           index_sta = find(possible_p>=x_range_sta(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))));
           index_50 = find(possible_p>=x_range_prop_50(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))));
           index_75 = find(possible_p>=x_range_prop_75(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))));
           index_100 = find(possible_p>=x_range_prop_100(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))));
           x_range_sta(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))  = possible_p(index_sta(1))+rand(1);
           x_range_prop_50(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))  = possible_p(index_50(1))+rand(1);
           x_range_prop_75(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))  = possible_p(index_75(1))+rand(1);
           x_range_prop_100(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))  = possible_p(index_100(1))+rand(1);
        end
    end
    for k = 1:size(EdgeTable.EndNodes,1)
       if traffic_path{s}(k+size(NodeTable.Name,1)) == 0
            x_range_prop_50(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
            x_range_prop_75(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
            x_range_prop_100(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
       else
            x_range_prop_50(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = x_range_prop_50(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))+rand(1);
            x_range_prop_75(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = x_range_prop_75(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))+rand(1);
            x_range_prop_100(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = x_range_prop_100(k+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))+rand(1);
       end
    end
end

f_sta{seed}(num_iteration,:) = ObjFunc_per_slice_dvfs_final(x_range_sta,traffic_path,coefficient_alpha,NodeTable,EdgeTable,num_slices);
f_prop_50{seed}(num_iteration,:) = ObjFunc_per_slice_dvfs_final(x_range_prop_50,traffic_path,coefficient_alpha,NodeTable,EdgeTable,num_slices);
f_prop_75{seed}(num_iteration,:)  = ObjFunc_per_slice_dvfs_final(x_range_prop_75,traffic_path,coefficient_alpha,NodeTable,EdgeTable,num_slices);
f_prop_100{seed}(num_iteration,:)  = ObjFunc_per_slice_dvfs_final(x_range_prop_100,traffic_path,coefficient_alpha,NodeTable,EdgeTable,num_slices);
f_opt{seed}(num_iteration,:)  = ObjFunc_per_slice_dvfs_final(x,traffic_path,coefficient_alpha,NodeTable,EdgeTable,num_slices);
%[T_sta, prob_sta] = delay_large_other_stra(x_range_min+1,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,replicas_destinations);
% [T_prop_50{seed}(num_iteration,:), prob_prop_50{seed}(num_iteration,:)] = delay_large_hyex(x_range_prop_50,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations);
% [T_prop_75{seed}(num_iteration,:), prob_prop_75{seed}(num_iteration,:)] = delay_large_hyex(x_range_prop_75,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations);
% [T_prop_100{seed}(num_iteration,:), prob_prop_100{seed}(num_iteration,:)] = delay_large_hyex(x_range_prop_100,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations);
% [T_opt{seed}(num_iteration,:), prob_opt{seed}(num_iteration,:)] = delay_large_hyex(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,replicas_destinations);
%prob_prop_50{seed}(num_iteration,:)= delay_large_hyex(x_range_prop_50,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations);
prob_prop_50{seed}(num_iteration,:) = delay_large_real_hyex(x_range_prop_50,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,routing_paths,prob_paths,replicas_destinations,num_components_slice);
prob_prop_75{seed}(num_iteration,:) = delay_large_real_hyex(x_range_prop_75,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,routing_paths,prob_paths,replicas_destinations,num_components_slice);
prob_prop_100{seed}(num_iteration,:) = delay_large_real_hyex(x_range_prop_100,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,routing_paths,prob_paths,replicas_destinations,num_components_slice);

%prob_prop_75{seed}(num_iteration,:) = delay_large_hyex(x_range_prop_75,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations);
%prob_prop_100{seed}(num_iteration,:) = delay_large_hyex(x_range_prop_100,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations);
prob_opt{seed}(num_iteration,:) = prob;


rep{seed} = replicas_destinations;
%save(['com_1s_50ms_5ms_slices_final_',num2str(seed),'_complexity_delay1e4_iter200_contolerance1-2_inter_',num2str(num_iteration),'_test_greedy_',num2str(greedy_num),'.mat'],'x','fval','res','f_sta','f_prop_50','f_prop_75','f_prop_100','f_opt','prob_opt','prob_prop_100','prob_prop_75','prob_prop_50','x_range_min','x_range_prop_50','x_range_prop_75','x_range_prop_100','prob_fail','rep','t_app');
save(['final_real_routing_1s_50ms_5ms_slices_final_Twaren.graphml_',num2str(seed),'_complexity_delay1e4_iter200_contolerance1-2_inter_',num2str(num_iteration),'_test_greedy_',num2str(greedy_num),'.mat'],'x','fval','res','f_sta','f_prop_50','f_prop_75','f_prop_100','f_opt','prob_opt','prob_prop_100','prob_prop_75','prob_prop_50','x_range_min','x_range_prop_50','x_range_prop_75','x_range_prop_100','prob_fail','rep','t_app','tEnd');

end
end
end

%objective function
function f = ObjFunc(x,traffic_path,coefficient_alpha,NodeTable,EdgeTable,num_slices)
f = 0;
%coefficient_nodes = 42.29;
coefficient_bw1 = 14.55/550;
coefficient_bw2 = 4.5;
coefficient_bw3 = 19.055;
%rho = zeros(num_slices,size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1));
count_rho = zeros(num_slices,size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1));
for i = 1:num_slices
    for j = 1:size(NodeTable.Name,1)
        if  x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) ~= 0
            count_rho(i,j) = coefficient_alpha(1,i)*traffic_path{i}(j);
            p_queue = 1/x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) *count_rho(i,j)*(5.2364*(x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/3e3)^2 - 13.4242 *x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/3e3 + 15.2484);
            f = f + p_queue;
            %rho(i,j) = count_rho(i,j)/x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)));
        end
    end
end



for i = 1:size(EdgeTable.EndNodes,1)
    temp2 = 0;
    for j = 1:num_slices
        temp2 = temp2 + x(i+size(NodeTable.Name,1)+(j-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)));
    end
    f = f + 0*(temp2==0) + (coefficient_bw2+coefficient_bw1*temp2).*(temp2>0&temp2<550)+(coefficient_bw3+0.0001*temp2).*(temp2>=550);
end
end

%%



