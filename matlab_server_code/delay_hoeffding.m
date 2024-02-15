function [ ceq,out,gradceq,gradc] = delay_hoeffding(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,replicas_destinations,x_range_max,t)
%arrival_app1 = random('Poisson', LambdaMicrophones,1,5);
%arrival_app2 = random('Poisson', LambdaMicrophones,1,5);

c = zeros(1,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices);
c_min = zeros(1,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices);
gradc = zeros((size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices,num_slices);
%t = zeros(1,num_slices);
for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)
        if  traffic_path{s}(i) == 0 || (ismember(i,replicas_destinations{s}{1})&&(~ismember(i,[replicas_destinations{s}{2:end}])))
            c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
            c_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (1/total_arrival(s) * traffic_path{s}(i)/(x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_alpha(i,s)-traffic_path{s}(i)))*1e4;
            c_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (1/total_arrival(s) * traffic_path{s}(i)/(x_range_max(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_alpha(i,s)-traffic_path{s}(i)))*1e4;
            if c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) < 0
                c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 1e10;
            end
            if c_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) < 0
                c_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 1e10;
            end
        end
    end
end
for s = 1:num_slices
    for i = 1:size(EdgeTable.EndNodes,1)
        if traffic_path{s}(i+size(NodeTable.Name,1)) == 0
            c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
            c_min(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (1/total_arrival(s) * traffic_path{s}(i+size(NodeTable.Name,1))/(x(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_beta(i,s)-traffic_path{s}(i+size(NodeTable.Name,1))))*1e4;
            c_min(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (1/total_arrival(s) * traffic_path{s}(i+size(NodeTable.Name,1))/(x_range_max(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_beta(i,s)-traffic_path{s}(i+size(NodeTable.Name,1))))*1e4;
            if c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) < 0
                c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 1e10;
            end
            if c_min(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) < 0
                c_min(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 1e10;
            end
        end
    end
end

%[t] = approxiamate_chebyshev(x,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations,x_range_max,step,round);

%[t,~,~] = approxiamate_t(dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations,x_range_max,step,round,upp_round,low_round);
%fprintf('t_app1: %f \n:', t(1));
%fprintf('t_app2: %f \n:', t(2));
%fprintf('t_app3: %f \n:', t(3));

for s = 1:num_slices
    out(s) =  sum(c(1+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1):(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*s)) + t(s) - dl(s);
end

ceq =[];
%if nargout > 2
    for s = 1:num_slices
        for i = 1:size(NodeTable.Name,1)
            if traffic_path{s}(i) == 0 || (ismember(i,replicas_destinations{s}{1})&&(~ismember(i,replicas_destinations{s}{2})))
                gradc(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)),s) = 0;
            else
                gradc(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)),s) = 1e4*traffic_path{s}(i)/(-total_arrival(s)*x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))^2/coefficient_alpha(i,s)+2*x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))*traffic_path{s}(i)*total_arrival(s)-coefficient_alpha(i,s)*total_arrival(s)*traffic_path{s}(i)^2);
            end
        end
    end
    for s = 1:num_slices
        for i = 1:size(EdgeTable.EndNodes,1)
            if traffic_path{s}(i+size(NodeTable.Name,1)) == 0
                 gradc(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)),s) = 0;
            else
                gradc(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)),s) = 1e4*traffic_path{s}(i+size(NodeTable.Name,1))/(-total_arrival(s)*x(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))^2/coefficient_beta(i,s)+2*x(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))*traffic_path{s}(i+size(NodeTable.Name,1))*total_arrival(s)-coefficient_beta(i,s)*total_arrival(s)*traffic_path{s}(i+size(NodeTable.Name,1))^2);
            end
        end
    end
    gradceq = [];
%end
end

%out2 = dl2 - sum(c(16:30));
%     out = [dl1 - sum(c(1:15)), dl2 - sum(c(16:30))];