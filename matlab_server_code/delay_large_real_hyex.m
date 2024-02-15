function prob = delay_large_real_hyex(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,routing_paths,prob_paths,replicas_destinations,num_components_slice)
%arrival_app1 = random('Poisson', LambdaMicrophones,1,5);
%arrival_app2 = random('Poisson', LambdaMicrophones,1,5);

c = zeros(1,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices);
flag = ones(1,num_slices);
%gradc = zeros((size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*num_slices,num_slices);
prob = zeros(1,num_slices);
for s = 1:num_slices
    for i = 1:size(NodeTable.Name,1)
        if  traffic_path{s}(i) == 0 || (ismember(i,replicas_destinations{s}{1})&&(~ismember(i,[replicas_destinations{s}{2:end}])))
            c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (1/total_arrival(s) * traffic_path{s}(i)/(x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_alpha(i,s)-traffic_path{s}(i)))/1e2;
            if c(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) < 0
                flag(s) = 0;
            end
        end
    end
end
for s = 1:num_slices
    for i = 1:size(EdgeTable.EndNodes,1)
        if traffic_path{s}(i+size(NodeTable.Name,1)) == 0
            c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = 0;
        else
            c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) = (1/total_arrival(s) * traffic_path{s}(i+size(NodeTable.Name,1))/(x(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/coefficient_beta(i,s)-traffic_path{s}(i+size(NodeTable.Name,1))))/1e2;
            if c(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) <0
                flag(s) = 0;
            end
        end
    end
end
q = {};
for s = 1:num_slices
    q{s} = [];
    for i = 1:size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)
        if c(i+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1)) ~= 0
            q{s}(i)= 1/c(i+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1));
        else 
            q{s}(i) = 0;
        end
        %out(s) =  sum(c(1+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1):(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*s)) - (dl(s)/(-log(prob_fail(s))));
    end
end
lam_path = {};
lam_route = {};
for s = 1:num_slices
    lam_route{s} = {};
    lam_path{s} = zeros(1, size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1));
    for i = 1:num_components_slice(s)-1 
        lam_route{s}{i} = {};
         start_point = replicas_destinations{s}{i};
         for start = 1:size(start_point,2)
             lam_route{s}{i}{start} = [];
             for j = 1:size(routing_paths{s}{i}{start},1)
                 %lam_route{s}{i}{start}(j) = [];
                 for k = 1:size(routing_paths{s}{i}{start},2)-1
                     val = find(EdgeTable.EndNodes(:,1) == routing_paths{s}{i}{start}(j,k) & EdgeTable.EndNodes(:,2) == routing_paths{s}{i}{start}(j,k+1));
                     lam_path{s}(val+size(NodeTable.Name,1)) = q{s}(val+size(NodeTable.Name,1));
                     lam_route{s}{i}{start}(j,k) =  lam_path{s}(val+size(NodeTable.Name,1));
                 end
             end
         end
    end
end
lam_total_route ={};
lam_total_prob = {};
%lam_final_route = {};
for s = 1:num_slices
    %5lam_final_route{s} = [];
    x_path = {};
    val = replicas_destinations{s};
    [x_path{1 : numel(val)}] = ndgrid(val{:});
    ret = reshape(cat(numel(val), x_path{:}), [], numel(val));
    lam_total_route{s} = {};
    for i = 1:size(ret,1)
        flag_exist = 0;
        lam_total_prob{s}(i) = 1;
        lam_total_route{s}{i} = [];
        for j = 2:size(ret,2)
            if ismember(q{s}(ret(i,j)),lam_total_route{s}{i}) == 0
                lam_total_route{s}{i} = [lam_total_route{s}{i} q{s}(ret(i,j))];
            end
        end
        %lam_final_route{s}(i,:) = lam_total_route{s}(i,:);
        for k = 1:size(ret,2)-1
            for m = 1:size(replicas_destinations{s}{k},2)
                for index = 1: size(routing_paths{s}{k}{m},1)
                    if ret(i,k) == routing_paths{s}{k}{m}(index,1) && ret(i,k+1) == routing_paths{s}{k}{m}(index,end)
                        %lam_final_route{s}(i,:) = [lam_final_route{s}(i,:) lam_route{s}{k}{m}(index,:)];
                        if ret(i,k) ~= ret(i,k+1)
                            lam_total_route{s}{i} = [lam_total_route{s}{i} lam_route{s}{k}{m}(index,:)];
                            lam_total_prob{s}(i) = lam_total_prob{s}(i)*prob_paths{s}{k}{m}(index);
                        end
                        flag_exist = flag_exist + 1;
                    end
                end
            end
        end
        if ismember(0,lam_total_route{s}{i})==1
            lam_total_route{s}{i} =[];
            lam_total_prob{s}(i) = 0;
        end
        if flag_exist ~= num_components_slice(s)-1
            lam_total_route{s}{i} =[];
            lam_total_prob{s}(i) = 0;
        end
    end
end

for s = 1:num_slices
    prob(s) = 1- hypo_multi_cdf(lam_total_route{s}, lam_total_prob{s},dl(s)/1e6);
end
% total_delay = zeros(1,num_slices);
% for s = 1:num_slices
%     for i = 1:length(q{s})
% %         if length(q{s}) == 1
% %             total_delay(s) = exp(-q{s}(i)*dl(s));
% %         else
%             delay_value = 1;
%             for j = 1:length(q{s})
%                 if j ~= i
%                     if abs(q{s}(j) - q{s}(i)) < 1
%                        q{s}(j) = q{s}(j)+1;
%                     end
%                     delay_value = delay_value * q{s}(j)/(q{s}(j)-q{s}(i));
%                     %delay_value
%                 end
%             end
%             total_delay(s) = total_delay(s) + delay_value*exp(-q{s}(i)*dl(s)/1e6);
%             %total_delay
% %   end
%     end
%     prob(s) = 1-total_delay(s);
%     if flag(s) == 0
%         prob(s) = 0;
%     end
% end 

% if sum( ( prob >= 0) & (prob <=1.0001)) ~= 3
% %%%%%cluster
% beta = {};
% number = {};
% %r = zeros(1,num_slices);
% slice_count = zeros(1,num_slices);
% for s = 1:num_slices
%     for i = 1:length(q{s})-1
%         for j = (i+1):length(q{s})
%             if abs(q{s}(j)-q{s}(i)) < 1 && (q{s}(j)~= q{s}(i))
%                 q{s}(j) = q{s}(i);
%             end
%         end
%     end
%     slice_count(s) = length(unique(q{s}));
%     r = unique(q{s});
%     for m = 1:slice_count(s)
%         beta{s}(m) = r(m);
%         number{s}(m) = length(find(q{s} == r(m)));
%     end
% end

% B = ones(1,num_slices);
% %R = zeros(1,num_slices);
% for s = 1: num_slices
%     for i = 1: slice_count(s)
%         B(s) = B(s)*beta{s}(i)^number{s}(i);
%     end
% end

% 
% prob_erlang = zeros(1,num_slices);
% for s = 1:num_slices
%     if prob(s) < 0 || prob(s) >= 1.0001 
%         syms lap
%         F = B(s)*1/lap;
%         for i = 1: slice_count(s)
%             F = F * (1/(lap+beta{s}(i))^number{s}(i));
%         end
%         f = ilaplace(F);
%         t = dl(s)/1e6;
%         prob_erlang(s) = eval(f);
%     else
%         prob_erlang(s) = prob(s);
%     end
% end
% 
% 
% prob = prob_erlang;
% 
% 
% psi = {};
% variable_i = {};
% for s = 1: num_slices 
%     if isnan(prob_erlang(s))
%         for k = 1: slice_count(s)
%             for l = 1: number{s}(k)
%                 omiga_2 = 0;
%                 for j = 0 : slice_count(s)
%                     variable_i{s}(j+1) = l - 1;
%                     if j ~= k
%                         omiga_2 = omiga_2 + variable_i{s}(j+1);
%                     end
%                     if j == 0
%                         sum_var= (dl(s)/1e6)^(-(1+variable_i{s}(j+1)));
%                     else
%                         sum_var= sum_var* factorial(variable_i{s}(j+1)+number{s}(j)-1)/(factorial(number{s}(j)-1)*factorial(variable_i{s}(j+1)))*((beta{s}(j)+dl(s)/1e6)^(-(number{s}(j)+variable_i{s}(j+1))));
%                     end
%                 end
%                 sum_var = sum_var *omiga_2;
%                 psi{s}{k}(l) = (-1)^(l-1)*factorial(l-1) * sum_var;
%                 R(s) = R(s) + (psi{s}{k}(l)*(-beta{s}(k))*(dl(s)/1e6)^(number{s}(k)-l)*exp(-beta{s}(k)*dl(s)/1e6))/(factorial(number{s}(k)-l)*factorial(l-1));
%             end
%         end
%         R(s) = 1- R(s)*B(s);
%         prob(s) = R(s);
%     end
% end
%prob = prob_erlang;
% for s = 1:num_slices
%     out(s) =  sum(c(1+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1):(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*s)) - (dl(s)/(-log(prob_fail(s))));
% end
end
