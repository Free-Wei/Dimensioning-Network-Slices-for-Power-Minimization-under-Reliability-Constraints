function paths = find_paths(num_components_slices,replicas_destinations,routing_paths,prob_paths,lam_path)
    i = num_components_slices;
    for n = 1 : size(replicas_destinations{s}{i},2)
        paths(n) = [find_paths(i,routing_paths,prob_paths,lam_path); [lam_path{s}{i}{n}]]; 
    end
    i = i-1;
end