function [G, EdgeTable, NodeTable]  = read_graph(file)
S = fileread(file);
nodes = regexp(S, 'node.*?id="(?<id>\d+)"', 'names');
edges = regexp(S, 'edge.*?source="(?<source>\d+)".*?target="(?<target>\d+)"', 'names');
all_ids = {nodes.id};
data_id_temp = str2num(char(all_ids))+1;
data_id = cellstr(string(data_id_temp'));
all_sources = {edges.source};
all_targets = {edges.target};
data_sources_temp = str2num(char(all_sources))+1;
data_sources = cellstr(string(data_sources_temp'));
data_target_temp = str2num(char(all_targets))+1;
data_targets = cellstr(string(data_target_temp'));
[source_found, s] = ismember(data_sources, data_id);
nfidx = find(~source_found);
if ~isempty(nfidx)
   error('Source ids not found in node list, starting from "%s"', edges(nfidx(1).source));
end
[target_found, t] = ismember(data_targets, data_id);
nfidx = find(~target_found);
if ~isempty(nfidx)
   error('Target ids not found in node list, starting from "%s"', edges(nfidx(1).target));
end
EdgeTable1 = table([s.', t.'], 'VariableNames', {'EndNodes'});
EdgeTable2 = table([t.', s.'], 'VariableNames', {'EndNodes'});
EdgeTable = unique([EdgeTable1; EdgeTable2],"rows","stable");
%EdgeTable = [EdgeTable1; EdgeTable2];
NodeTable = table(data_id.', 'VariableNames',{'Name'});
G = digraph(EdgeTable,NodeTable);
plot(G)
end