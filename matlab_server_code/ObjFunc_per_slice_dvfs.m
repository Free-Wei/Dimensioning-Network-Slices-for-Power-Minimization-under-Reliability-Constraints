function f_final = ObjFunc_per_slice_dvfs(x,traffic_path,coefficient_alpha,NodeTable,EdgeTable,num_slices)
f = zeros(num_slices,(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)));
%coefficient_nodes = 42.29;
coefficient_bw1 = 14.55/550;
coefficient_bw2 = 4.5;
coefficient_bw3 = 19.055;

node_allocation = zeros(1,size(NodeTable.Name,1));
rho = zeros(num_slices,size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1));
count_rho = zeros(num_slices,size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1));
for i = 1:num_slices
    for j = 1:size(NodeTable.Name,1)
        if  x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1))) ~= 0
            count_rho(i,j) = coefficient_alpha(1,i)*traffic_path{i}(j);
            rho(i,j) = count_rho(i,j)/x(j+(i-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)));
        end
    end
end

for i = 1:size(NodeTable.Name,1)
    temp = 1;
    f_node = 0;
    for s = 1:num_slices
        node_allocation(i) = node_allocation(i) + x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)));
        temp = temp * (1-rho(s,i));
    end
    f_node = 0*(node_allocation(i)==0) + (0.956^2*temp*12)*((node_allocation(i)/NodeTable.Cpu(i)<= 0.375)&(node_allocation(i)/NodeTable.Cpu(i) > 0)) +...
        (1.036^2*temp*16)*((node_allocation(i)/NodeTable.Cpu(i)<= 0.5)& (node_allocation(i)/NodeTable.Cpu(i)> 0.375)) +...
        (1.164^2*temp*20)*((node_allocation(i)/NodeTable.Cpu(i)<= 0.625)& (node_allocation(i)/NodeTable.Cpu(i)> 0.5)) +...
        (1.276^2*temp*24)*((node_allocation(i)/NodeTable.Cpu(i)<= 0.75)& (node_allocation(i)/NodeTable.Cpu(i)> 0.625)) +...
        (1.42^2*temp*28)*((node_allocation(i)/NodeTable.Cpu(i)<= 0.875)& (node_allocation(i)/NodeTable.Cpu(i)> 0.75)) +...
        (1.484^2*temp*32)*((node_allocation(i)/NodeTable.Cpu(i)<= 1)& (node_allocation(i)/NodeTable.Cpu(i)> 0.875)) ;
    for s = 1:num_slices
        if node_allocation(i) == 0
            f(s,i) = 0;
        else
            f(s,i) = x(i+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/node_allocation(i) * f_node;
        end
    end
end

for i = 1:size(EdgeTable.EndNodes,1)
    temp2 = 0;
    for j = 1:num_slices
        temp2 = temp2 + x(i+size(NodeTable.Name,1)+(j-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)));
    end
    f_total_per_link = 0*(temp2==0) + (coefficient_bw2+coefficient_bw1*temp2).*(temp2>0&temp2<550)+(coefficient_bw3+0.0001*temp2).*(temp2>=550);
    for s = 1:num_slices
        if temp2 == 0
            f(s,i+size(NodeTable.Name,1)) = 0;
        else
            f(s,i+size(NodeTable.Name,1)) = x(i+size(NodeTable.Name,1)+(s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)))/temp2 * f_total_per_link;
        end
    end
end
f_final = sum(f,2);
end