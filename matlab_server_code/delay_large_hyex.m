function prob = delay_large_hyex(x,dl,traffic_path,coefficient_alpha,coefficient_beta,NodeTable,EdgeTable,total_arrival,num_slices,~,replicas_destinations)
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
    count = 0;
    for i = 1:size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)
        if c(i+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1)) ~= 0
            count = count + 1;
            q{s}(count)= 1/c(i+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1));
        end
        %out(s) =  sum(c(1+(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*(s-1):(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))*s)) - (dl(s)/(-log(prob_fail(s))));
    end
end
total_delay = zeros(1,num_slices);
for s = 1:num_slices
    for i = 1:length(q{s})
%         if length(q{s}) == 1
%             total_delay(s) = exp(-q{s}(i)*dl(s));
%         else
            delay_value = 1;
            for j = 1:length(q{s})
                if j ~= i
                    if abs(q{s}(j) - q{s}(i)) < 1
                       q{s}(j) = q{s}(j)+1;
                    end
                    delay_value = delay_value * q{s}(j)/(q{s}(j)-q{s}(i));
                    %delay_value
                end
            end
            total_delay(s) = total_delay(s) + delay_value*exp(-q{s}(i)*dl(s)/1e6);
            %total_delay
%   end
    end
    prob(s) = 1-total_delay(s);
    if flag(s) == 0
        prob(s) = 0;
    end
end 

if sum( ( prob >= 0) & (prob <=1.0001)) ~= 3
%%%%%cluster
beta = {};
number = {};
%r = zeros(1,num_slices);
slice_count = zeros(1,num_slices);
for s = 1:num_slices
    for i = 1:length(q{s})-1
        for j = (i+1):length(q{s})
            if abs(q{s}(j)-q{s}(i)) < 1 && (q{s}(j)~= q{s}(i))
                q{s}(j) = q{s}(i);
            end
        end
    end
    slice_count(s) = length(unique(q{s}));
    r = unique(q{s});
    for m = 1:slice_count(s)
        beta{s}(m) = r(m);
        number{s}(m) = length(find(q{s} == r(m)));
    end
end

B = ones(1,num_slices);
%R = zeros(1,num_slices);
for s = 1: num_slices
    for i = 1: slice_count(s)
        B(s) = B(s)*beta{s}(i)^number{s}(i);
    end
end


prob_erlang = zeros(1,num_slices);
for s = 1:num_slices
    if prob(s) < 0 || prob(s) >= 1.0001 
        syms lap
        F = B(s)*1/lap;
        for i = 1: slice_count(s)
            F = F * (1/(lap+beta{s}(i))^number{s}(i));
        end
        f = ilaplace(F);
        t = dl(s)/1e6;
        prob_erlang(s) = eval(f);
    else
        prob_erlang(s) = prob(s);
    end
end


prob = prob_erlang;
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
end
% if nargout > 2
%     for s = 1:num_slices
%         for i = 1:size(NodeTable.Name,1)
%             if traffic_path{s}(i) == 0
%                 gradc(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)),s) = 0;
%             else
%                 gradc(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)),s) = traffic_path{s}(i)/(-total_arrival(s)*x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))^2/coefficient_alpha(i,s)+2*x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))*traffic_path{s}(i)*total_arrival(s)-coefficient_alpha(i,s)*total_arrival(s)*traffic_path{s}(i)^2);
%             end
%         end
%     end
%     for s = 1:num_slices
%         for i = 1:size(EdgeTable.EndNodes,1)
%             if traffic_path{s}(i+size(NodeTable.Name,1)) == 0
%                  gradc(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)),s) = 0;
%             else
%                 gradc(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)),s) = traffic_path{s}(i+size(NodeTable.Name,1))/(-total_arrival(s)*x(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))^2/coefficient_beta(i,s)+2*x(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))*traffic_path{s}(i+size(NodeTable.Name,1))*total_arrival(s)-coefficient_beta(i,s)*total_arrival(s)*traffic_path{s}(i+size(NodeTable.Name,1))^2);
%             end
%         end
%     end
%     gradceq = [];
% end
% end

%out2 = dl2 - sum(c(16:30));
%     out = [dl1 - sum(c(1:15)), dl2 - sum(c(16:30))];