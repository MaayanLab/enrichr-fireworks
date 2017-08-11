          
clear; clc; close;
files=repmat(string(), [1 11]);
% file size contains the length of each library in the order that it is
% read in files, this is used to assigning theright library to each gene
% set
filesize=[84;96;2192;2192;2918;2918;967;93;207;30;315]; 

files(1,1)='Human_Gene_Atlas.txt'
files(1,2)='Mouse_Gene_Atlas.txt'
files(1,3)='Allen_Brain_Atlas_up.txt'
files(1,4)='Allen_Brain_Atlas_down.txt'
files(1,5)='GTEx_Tissue_Sample_Gene_Expression_Profiles_up.txt'
files(1,6)='GTEx_Tissue_Sample_Gene_Expression_Profiles_down.txt'
files(1,7)='Cancer_Cell_Line_Encyclopedia.txt'
files(1,8)='NCI-60_Cancer_Cell_Lines.txt'
files(1,9)='Tissue_Protein_Expression_from_ProteomicsDB.txt'
files(1,10)='Tissue_Protein_Expression_from_Human_Proteome_Map.txt'
files(1,11)='ESCAPE.txt' 
%             
            fout=fopen('alltext.txt','w');
           for i=1:length(files)
               ID=fopen(files(1,i),'r');
               temp=fread(ID,'uint8');
               fwrite(fout,temp,'uint8');
               fclose(ID);
               
           end
           ID=fopen('alltext.txt','r');
           S=fscanf(ID,'%c');
           fclose(ID);
            fclose('all');
           % read in human and mouse gene sets
           %compare all gmt file genes to these sets
           %if gene exists in one of the two sets it can be kept
          %this is to reduce non-genes from gmt typos
           idmouse=fopen('mouse_gene_symbols.txt','r');
           mouse=fscanf(idmouse,'%c');
           mouse=strsplit(mouse,{'\n','\r','\v'});
           mouse=upper(mouse);
           fclose(idmouse);
           idhuman=fopen('human_gene_symbols.txt','r');
           human=fscanf(idhuman,'%c');
           human=strsplit(human,{'\n','\r','\v'});
           human=upper(human);
           fclose(idhuman);           
            maxnumgenes=25000;
            %split string based on new line
            c=strsplit(S,{'\n','\r','\v'});
            c=transpose(c);
            names=repmat(string(),[length(c) 1]);
            genes=repmat(string(),[length(c) maxnumgenes]);
            allgenes=[];
            genelength=zeros(length(c),1);
         
%take the gene set rows, for each row split into genes, standardize
%uppercase, make sure it belongs in either human or mouse set, add new genes
%to allgenes set and add corresponding info for name matrix and gene matrix

            for i=1:length(c)
                %split gene set based on space
                B=strsplit(string(c(i)));
                B_genesonly=B(2:end);
                B_genesonly=upper(B_genesonly);
                B_genesonly=strtok(B_genesonly,',');
                % convert this to a cell array, if that array is not empty,
                % delete the whole element
                B_genesonly(ismember(B_genesonly,[mouse human])==0)=[];
                names(i,1)=B(1);
                genes(i,1:length(B_genesonly))=B_genesonly;
                genelength(i)=length(B_genesonly);
                subset={allgenes,B_genesonly};
                allgenes=unique(cat(2,subset{:}));
            end

%%
% create the name matrix by giving each gene set the correctlibrary number
 names(:,2)=zeros(length(names),1);
counter=1;
for i=1:length(files)
    for j=1:filesize(i)
         names(counter,2)=i;
          counter=counter+1;        
    end
end
namescopy=names;
genescopy=genes;

%% get rid of gene sets under a certain length
% genethresh=prctile(genelength,25,1);
% [r,c]=find(genelength<genethresh);
% genelength(r,:)=[];
% genes(r,:)=[];
% names(r,:)=[];

     %% 
     %helpgenes is a matrix where each row is a gene set, columns are the
     %set of all genes, if a gene is present that element has a 1
            len=size(genes);
            len=len(1);            
            helpgenes=zeros(len,length(allgenes));

            for i=1:len
                %get number of actual genes, not spaces
                for j=1:genelength(i)
                    ind=find(allgenes==string(genes(i,j)));
                    helpgenes(i,ind)=1;
                end
            end
%             
            

%%  
%send help genes to calculate the jaccard to get an adj mat between all
%gene sets
         find_intersect=squareform(pdist(helpgenes,'jaccard'));

%%
%pdist gives the inverse values, use 1-find_intersect to get actual adj mat
  [r,c]=find(isnan(find_intersect));
 find_intersect(r,c)=1;
 D=1-find_intersect;
 D=D-eye(length(D));
 
 %%
  %% save adjmat to txt file
% completedpathways=sparse(D);
%  [i,j,val] = find(completedpathways);
% data_dump3 = [i,j,val];
%  completedpathways = spconvert( data_dump3 );
%  save -ascii DiseaseJun31adjmat.txt data_dump3
%  fid = fopen('DiseaseJun31adjmat.txt','w');
% fprintf( fid,'%d %d %f\n', transpose(data_dump3) );
% fclose(fid);
%% if you have previously calculated adjmat, you can open it here

% %%take txt files to matrices
% textname='CellTypeJun31adjmat.txt';
% id2=fopen(textname,'r');
% path=fscanf(id2,'%f',[3 Inf]);
% fclose(id2);
% path=path';
% PATH=spconvert(path);
% PATH=full(PATH);
%%
%PATH is the basic name used for the adjacency matrix 
PATH=D;
%% if you would like to get rid of gene sets that have below a certain number of genes
% genethresh=prctile(genelength,25,1);
% [r,c]=find(genelength<genethresh);
% genelength(r,:)=[];
% genes(r,:)=[];
% names(r,:)=[];
% PATH(r,:)=[];
% PATH(:,r)=[];]]
%% if you would like to delete links that are below a certain weight
%threshold=prctile(reshape(PATH,1,[]),99);
%PATH(PATH<threshold)=0;

%% KNN, keep the top three highest links for every node to reduce density
PATH=D;

for i=1:length(PATH)
    a=sort(PATH(i,:),'descend');
    p=PATH(i,:);
    p(p<a(3))=0;
    PATH(i,:)=p;
end
for i=1:length(PATH)
    for j=i:length(PATH)
    if PATH(i,j)~=PATH(j,i)
        PATH(i,j)=max([PATH(i,j) PATH(j,i)]);
        PATH(j,i)=PATH(i,j);
    end
    end
end
%% delete unconnected components 
%first convert graph into unweighted in order to use connected components
%function
%connectedcomponents is a vector that assigns each node to its connected
%group
%use hist function to tally the number of nodes in each group, delete
%groups that have 10 nodes assigned to them or less, assigned in n
%delete those nodes from the PATH graph, and from the names and genes
%matrices
unweighted=PATH;
unweighted(unweighted>0)=1;
unweighted_graph=graph(unweighted);
connectedcomponents=conncomp(unweighted_graph);
group_tally=hist(connectedcomponents,max(connectedcomponents));
min_groups=find(group_tally<10);
n=zeros(length(PATH),1);
nindex=1;
for i=1:length(min_groups)
r=min_groups(i);
nodes_for_group_r=find(connectedcomponents==r);
for j=1:length(nodes_for_group_r)
%append node onto end of vector
n(nindex)=nodes_for_group_r(j);
nindex=nindex+1;
end
end

n(n==0)=[];
PATH(n,:)=[];
PATH(:,n)=[];
names(n,:)=[];
genes(n,:)=[];

unweighted=PATH;
unweighted(unweighted>0)=1;
unweighted_graph=graph(unweighted);
%plot the degree distribution of nodes using the updated unweighted graph
connections=sum(unweighted);
tot=connections;
cumu_connections=zeros(1,length(PATH));
 for i=1:length(PATH)-1
    x=tot(i)+1; 
    if (x>0)
    cumu_connections(x)=cumu_connections(x)+1 ;
    end
 end  
scatter(log(1:max(connections)),log(cumu_connections(1:max(connections))))      
xlabel('degree of node')
ylabel('number of nodes')

% save files as a matlab file in order to apply python layout, and as gml
% file to appy cytoscape layout

%save DisKNN.mat PATH -v7.3
%matToGML(PATH,'fileName','DisKNN');

%%
%add index to names and genes matrices and write to text files
names2=names;
names3=strrep(names2,',','_');
numberofgenes=size(names);
numberofgenes=numberofgenes(1);
nameindex=transpose(string(1:numberofgenes));
nameindex(:,2:3)=names3;
fid=fopen('CellTypeKNNnames.txt','w');
 fprintf(fid, '%s,%s,%s\r\n',transpose(nameindex));

genesfill2=fillmissing(genes,'constant',string(string('')));
%convert row of string elements into one string element
combinedgenes=repmat(string(),[numberofgenes 1]);
for i=1:numberofgenes
combinedgenes(i)=strjoin(genesfill2(i,:));
end
%add index
index=string(transpose([1:numberofgenes]));
combinedgenes2=[index,combinedgenes];
%write to proper file
fid=fopen('CellTypeKNNgenes.txt','w');
fprintf(fid,'%s %s\r\n',transpose(combinedgenes2));
