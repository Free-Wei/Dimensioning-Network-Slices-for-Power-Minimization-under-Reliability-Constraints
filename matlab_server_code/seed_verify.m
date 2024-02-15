function out = seed_verify(seed)
flag_cycle = 1;
while flag_cycle == 1
flag_cycle = 0;
rng(seed);
[G, EdgeTable, NodeTable] = read_graph("Twaren.graphml");
Cpu = ones(size(NodeTable.Name,1),1)*1.032e4; 
% 10GHz
Bandwidth = ones(size(EdgeTable.EndNodes,1),1)*1e4;  
%Bandwidth = ((b-a).*rand(size(EdgeTable.EndNodes,1),1)+a)*1e4;
NodeTable = addvars(NodeTable,Cpu);
EdgeTable = addvars(EdgeTable,Bandwidth);
num_slices = 3;
min_num_components = 2;
max_num_components = 5;
max_number_of_replicas = 2;
max_number_of_ingress = 2;
num_components_slice = [randi([min_num_components,max_num_components],1,num_slices-1),randi([2,3],1,1)];
%dl = [1,0.05,0.005]*1e4;
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
arrival_rate = [5000 100 3000]+ [1000 80 200]*5; 
%p1 = [0.5,0.5];
p_route = 0.5;
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
                if (isempty(shortestpath(G, start_point(start), finish_point(finish)))==1 && start_point(start)~=finish_point(finish))
                    flag_cycle = 1;
                end
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
if flag_cycle == 1
    seed = seed+1;
end
end
out = seed;
end