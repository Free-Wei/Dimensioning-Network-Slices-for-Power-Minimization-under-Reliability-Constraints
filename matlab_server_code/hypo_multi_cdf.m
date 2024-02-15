function res = hypo_multi_cdf(lam_total_route, lam_total_prob,T)
out = [];
for k = 1:length(lam_total_route)
    if isempty(lam_total_route{k}) == 0
    out_one_route = [];
    for i = 1:length(lam_total_route{k})
        temp = 1;
        for j = 1 : length(lam_total_route{k})
            if j ~= i
                if abs(lam_total_route{k}(j) - lam_total_route{k}(i)) < 1
                       lam_total_route{k}(j) = lam_total_route{k}(j)+1;
                end
                temp = lam_total_route{k}(j)/(lam_total_route{k}(j)-lam_total_route{k}(i)) *temp; 
            end
        end
        out_one_route(i) = temp *exp(-lam_total_route{k}(i)*T);
    end
    out(k) = lam_total_prob(k) * sum(out_one_route);
    end
end
res = sum(out);
end