function [out,gradc] = ineq_fn(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations,x_range_max,round)
%arrival_app1 = random('Poisson', LambdaMicrophones,1,5);
%arrival_app2 = random('Poisson', LambdaMicrophones,1,5);

c = zeros(1,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices);
c_min = zeros(1,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices);
gradc = zeros((size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices,num_slices);
t = zeros(1,num_slices);
for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)
        if  traffic_path{s}(i) == 0 || (ismember(i,replicas_destinations{s}{1})&&(~ismember(i,replicas_destinations{s}{2})))
            c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
            c_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (1/total_arrival(s) * traffic_path{s}(i)/(x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_alpha(i,s)-traffic_path{s}(i)))*1e4;
            c_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (1/total_arrival(s) * traffic_path{s}(i)/(x_range_max(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_alpha(i,s)-traffic_path{s}(i)))*1e4;
            if c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) < 0
                c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = -1e9;
                c_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = -1e9;
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
        end
    end
end
q_upper = {};
q_lower = {};
rate = {};
for s = 1:num_slices
    count = 0;
    for i = 1:size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)
        if c(i+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1)) ~= 0
            count = count + 1;
            q_upper{s}(count)= c(i+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1));
            rate{s}(count) =  traffic_path{s}(i);
            q_lower{s}(count)= c_min(i+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1));
        end
        %out(s) =  sum(c(1+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1):(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*s)) - (dl(s)/(-log(prob_fail(s))));
    end
end
upper_bound = {};
lower_bound = {};
result_hoeffding = zeros(1,num_slices);
for s = 1:num_slices
    total_rate = sum(rate{s});
    for i = 1:length(q_upper{s})
            %upper_bound{s}(i) = q_upper{s}(i) + 0.5*q_upper{s}(i)^2;
            upper_bound{s}(i) = dl(s)*rate{s}(i)/total_rate/round(s);
            lower_bound{s}(i) = q_lower{s}(i);
            result_hoeffding(s) = result_hoeffding(s) + (upper_bound{s}(i)-lower_bound{s}(i))^2;
        %delay_value = 1;
        %for j = 1:length(q{s})
            %%if j ~= i
            %    delay_value = delay_value * (q{s}(j)*exp(-q{s}(i)*dl(s)))/(q{s}(j)-q{s}(i));
            %end
        %end
        %total_delay = total_delay + delay_value;
    end
    if length(q_upper{s}) == 1
        t(s) = -log(prob_fail(s))*1/q_upper{s}(1);
    else
    %out(s) = -total_delay+prob_fail(s);
        t(s) = sqrt(-log(prob_fail(s))*result_hoeffding(s)/2);
    end
end 

for s = 1:num_slices
    out(s) =  sum(c(1+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1):(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*s)) - (dl(s)-t(s));
end

%ceq =[];
if nargout > 2
    for s = 1:num_slices
        for i = 1:size(NodeTable.Name,1)
            if traffic_path{s}(i) == 0
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
    %gradceq = [];
end
%inequlity capacity
out1 = zeros(1,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices);
out2 = zeros(1,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices);
for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)
        out1 = x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) - NodeTable.Cpu(i);

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
        end
    end
end
end

%out2 = dl2 - sum(c(16:30));
%     out = [dl1 - sum(c(1:15)), dl2 - sum(c(16:30))];