clear
tic
%% Inputs

V=1:3;
E=[2 1 1 2 2 2 3 2 2 3 1 3];
Wadj=[0 3 2 0 0;2 0 1 0 0;0 1 0 0 0;0 0 0 0 0;0 0 0 0 0];
ntest=0; 
% Write "0" for "ntest" if minimal expansion is intended; otherwise, enter an array 
% with the same # elements as "V" to check for feasibility or to generate non-minimal network.

%% Conversion to script's reading structure

s=diag(Wadj)'; % Create 's' array.
c=E;
for i=1:2:length(E)-1 % Remove self-loops from 'c' array.
    if c(i)==c(i+1)
        c(i:i+1)=0;
    end
end
c=nonzeros(c)';
for i=1:2:length(c) % Remove reverse edge pairs from 'c' array.
    if i+3<=length(c)&&all(c(i:i+1)==c(i+3:-1:i+2))
        c(i+2:i+3)=0;
    end
end
c=nonzeros(c)'; % Create 'c' array.
f=zeros(1,length(c));
for i=1:2:length(c)
    f(i)=Wadj(c(i+1),c(i));
    f(i+1)=Wadj(c(i),c(i+1)); % Create 'f' array.
end

% Notes: s = number of self loops per quotient node, f = vector pairs of edge weights.
% c = quotient node pairs with respect to array f. The index of c is the node that generates a given number of edges given by the same index in f. 
% These edges travel to the other quotient node in that pair of selected nodes, and vice versa.

%% Step 1  - Number of expanded nodes per quotient node (factoring in both non self and self-loops)

% - j and j1 are dummy indeces
% - ind is an indexing array

L=length(f);

n=zeros(1,max(c)); % Can be max(c) OR length(s)
for i=1:max(c) % Creates the n vector: The number of expanded nodes per quotient node.
    ind=ones(1,L);
    for j=1:L
        if c(j)==i
            ind(j)=f(j);
        end
    end
    if s(i)>0
        ind(j+1)=s(i)+1;
    end
    n(i)=max(ind);
end

%% Inputs feasibility

if ntest~=0
    test=ntest<n;
    if sum(test)>0    % If feasibility script is commented out, error message stems from "randperm" functions, indicating the graph is infeasible.
    disp('Graph is infeasible')
    return
    end
    n=ntest;
end

%% Step 2 - Assign the non-self-loop edges (random assignment)

% - pind and cind are the indexes for the parent and child nodes,
% respectively.
% - parent and child are the number of expanded parent and child nodes,
% respectively.
% - A is the adjacency matrix

A=zeros(sum(n));

for i=1:max(c) % Assigns the directed edges in the entire adjacency matrix, A.
    ind=zeros(1,L);
    for j=1:2:L
        pind=0;
        cind=0;
        if c(j)==i % Selects index of parent and child node in pair. If neither index is the target node, skips this.
            pind=j;
            cind=j+1;
            ind(j:j+1)=f(j:j+1);
        elseif c(j+1)==i
            pind=j+1;
            cind=j;
            ind(j:j+1)=f(j:j+1);
        end
        if pind~=0 && ind(pind)~=0 % If the parent node creates no edges to child or the target node isn't contained (based on above), skips this.
            child=n(c(cind));
            parent=n(c(pind));
            paredge=f(pind);
            padj=sum(n(1:c(pind)-1))+1:sum(n(1:c(pind))); % Indices for parent nodes in adjacency matrix.
            cadj=sum(n(1:c(cind)-1))+1:sum(n(1:c(cind))); % Indices for child nodes in adjacency matrix.
            if sum(sum(A(cadj,padj)>0))>0 % Checks for the occurrence of a repeated edge.
                disp('Error, repeated edge detected')
            end
            for j1=1:length(cadj) % j1 is the selected index of cadj
                A(cadj(j1),randperm(length(padj),f(pind))+min(padj)-1)=1; % For each child index, randomly assign unrepeated combination of parents with number of edges needed.
            end

        end
    end
end

%% Step 3 - Assign the self-loop edges

for i=1:max(c)
    if s(i)>0
        sadj=sum(n(1:i-1))+1:sum(n(1:i)); % sadj is the array of self-loop orbit nodes.
        for j=1:n(i)
            msadj=sadj; %msadj is the modifiable dummy array of sadj.
            msadj(msadj==msadj(j))=[];
            psadj=msadj(randperm(length(msadj),s(i)));
            A(sadj(j),psadj)=1;
        end
    end
end

%% Step 4 - Graph it and report in-degree

a=plot(flipedge(digraph(A)),MarkerSize=8);
color=hsv(length(n));
nodePartition(1:sum(n),1)=0;
for i=1:length(n)
    partitions=sum(n(1:i-1))+1:sum(n(1:i));
    highlight(a,partitions,'NodeColor',color(i,:))
    nodePartition(partitions,1)=i;
end
in_degree=sum(A,2);
Gephi_Colors=255*color;

%% Step 5 - Save node and edge lists to CSV

NNumber=(1:width(A));
zeroedgeID=zeros([length(nonzeros(A)) 2]);
edgeList={zeroedgeID};
k=0;
for i = 1:size(A,2)
    for j = 1:size(A,1)
        if A(i,j) ~= 0
            k=k+1;
            edgeList(k,[1 2])={NNumber(j) NNumber(i)};
        end
    end
end
nodeID=NNumber';
nodeList=num2cell([nodeID nodeID nodePartition]);

E=cell2table(edgeList); % Convert to table
N=cell2table(nodeList);
E.Properties.VariableNames=["Source","Target"];
N.Properties.VariableNames=["ID","Label","Partition"];
writetable(E, 'edge_list.csv'); % Write to CSV
writetable(N, 'node_list.csv');
disp('Node and list saved to node_list.csv and edge_list.csv, respectively.');
toc
