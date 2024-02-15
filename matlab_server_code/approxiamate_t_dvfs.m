function [t,round1,step] = approxiamate_t_dvfs(dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,prob_fail,replicas_destinations,x_range_max,step,round)

c_min = zeros(1,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices);
%upp_round1 = upp_round;
%low_round1 = low_round;
round1 = round;
for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)
        if  traffic_path{s}(i) == 0 || (ismember(i,replicas_destinations{s}{1})&&(~ismember(i,[replicas_destinations{s}{2:end}])))
            c_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            c_min(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (1/total_arrival(s) * traffic_path{s}(i)/(x_range_max(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_alpha(i,s)-traffic_path{s}(i)))*1e4;
        end
    end
end
for s = 1:num_slices
    for i = 1:size(EdgeTable.EndNodes,1)
        if traffic_path{s}(i+size(NodeTable.Name,1)) == 0
            c_min(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            c_min(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (1/total_arrival(s) * traffic_path{s}(i+size(NodeTable.Name,1))/(x_range_max(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_beta(i,s)-traffic_path{s}(i+size(NodeTable.Name,1))))*1e4;
        end
    end
end
q_lower = {};
%rate = {};
count = zeros(1,num_slices);
for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)
        if c_min(i+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1)) ~= 0
            count(s) = count(s) + 1;
            %rate{s}(count(s)) =  traffic_path{s}(i);
            q_lower{s}(count(s))= c_min(i+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1));
        end
        %out(s) =  sum(c(1+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1):(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*s)) - (dl(s)/(-log(prob_fail(s))));
    end
end
upper_bound = {};
lower_bound = {};
result_hoeffding = zeros(1,num_slices);
t = zeros(1,num_slices);
for s = 1:num_slices
    %total_rate = sum(rate{s});
    for i = 1:count(s)
        upper_bound{s}(i) = dl(s) * sqrt(-32  /(25 * count(s) * log(prob_fail(s))))/round1(s);
        lower_bound{s}(i) = 0;
        result_hoeffding(s) = result_hoeffding(s) + (upper_bound{s}(i)-lower_bound{s}(i))^2;
    end
    %round1(s) = (upp_round1(s)-low_round1(s))/2+low_round1(s); 
    t(s) = sqrt(-log(prob_fail(s))*result_hoeffding(s)/2);
end
%fprintf('round_next: %f \n:', round1);
fprintf('estimated_t: %f \n:', t);
% for s = 1:num_slices
%     while t(s) > dl(s)
%         result_hoeffding(s) = 0;
%         low_round1(s) = round1(s);
%         round1(s) = (upp_round1(s)-low_round1(s))/2 + low_round1(s);
%         %total_rate = sum(rate{s});
%         for i = 1:count(s)
%             %upper_bound{s}(i) = q_upper{s}(i) + 0.5*q_upper{s}(i)^2;
%             upper_bound{s}(i) = 32 * dl(s) /(25 * count(s) * log(prob_fail(s)))/round1(s);
%             %                 if upper_bound{s}(i) >= -q_upper{s}(i)*log(prob_fail(s))
%             %                     upper_bound{s}(i) = -q_upper{s}(i)*log(prob_fail(s));
%             %                 end
%             lower_bound{s}(i) = q_lower{s}(i);
%             result_hoeffding(s) = result_hoeffding(s) + (upper_bound{s}(i)-lower_bound{s}(i))^2;
%             %delay_value = 1;
%             %for j = 1:length(q{s})
%             %%if j ~= i
%             %    delay_value = delay_value * (q{s}(j)*exp(-q{s}(i)*dl(s)))/(q{s}(j)-q{s}(i));
%             %end
%             %end
%             %total_delay = total_delay + delay_value;
%         end
%         %             if length(q_upper{s}) == 1
%         %                 t(s) = -log(prob_fail(s))*1/q_upper{s}(1);
%         %             else
%         %out(s) = -total_delay+prob_fail(s);
%         t(s) = sqrt(-log(prob_fail(s))*result_hoeffding(s)/2);
%     end
% end

end