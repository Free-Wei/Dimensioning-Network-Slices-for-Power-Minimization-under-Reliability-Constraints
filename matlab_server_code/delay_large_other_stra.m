function [T_cost,prob_succ] = delay_large_other_stra(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,replicas_destinations)
%arrival_app1 = random('Poisson', LambdaMicrophones,1,5);
%arrival_app2 = random('Poisson', LambdaMicrophones,1,5);

c = zeros(1,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices);
T_cost = zeros(1,num_slices);
prob_succ = zeros(1,num_slices);
for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)
        if traffic_path{s}(i) == 0 || (ismember(i,replicas_destinations{s}{1})&&(~ismember(i,replicas_destinations{s}{2})))
            c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 1/total_arrival(s) * traffic_path{s}(i)/(x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_alpha(i,s)-traffic_path{s}(i));
            if c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) < 0
                c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = -1e9;
            end
        end
    end
end
for s = 1:num_slices
    for i = 1:size(EdgeTable.EndNodes,1)
        if traffic_path{s}(i+size(NodeTable.Name,1)) == 0
            c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 1/total_arrival(s) * traffic_path{s}(i+size(NodeTable.Name,1))/(x(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_beta(i,s)-traffic_path{s}(i+size(NodeTable.Name,1)));
            if c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) < 0
                c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = -1e9;
            end
        end
    end
end
for s = 1:num_slices
    T_cost(s) =  sum(c(1+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1):(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*s));
    if 1-exp(-dl(s)/T_cost(s)) < 0
        prob_succ(s) = 0;
    else
        prob_succ(s) = 1-exp(-dl(s)/T_cost(s));
    end
end
end