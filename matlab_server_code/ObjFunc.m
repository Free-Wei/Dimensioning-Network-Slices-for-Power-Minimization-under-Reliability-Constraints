function f = ObjFunc(x,NodeTable,EdgeTable,num_slices)
f = 0;
coefficient_nodes = 42.29;
coefficient_bw1 = 14.55/550;
coefficient_bw2 = 4.5;
coefficient_bw3 = 19.055;

for s = 1:num_slices
    temp = (s-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1));
    for i = 1:size(NodeTable.Name,1)
        f = f + coefficient_nodes/1e6*(x(i+temp));
    end
end
for i = 1:size(EdgeTable.EndNodes,1)
    temp2 = 0;
    for j = 1:num_slices
        temp2 = temp2 + x(i+size(NodeTable.Name,1)+(j-1)*(size(NodeTable.Name,1)+size(EdgeTable.EndNodes,1)));
    end
    f = f +(0)*(temp2==0) + (coefficient_bw2+coefficient_bw1*temp2).*(temp2>0&temp2<550)+(coefficient_bw3+0.0001*temp2).*(temp2>=550);
end
end
