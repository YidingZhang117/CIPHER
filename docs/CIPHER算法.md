# CIPHER算法

## 数据准备

1. 基因网络构建：ppi.txt 中描述了基因对应蛋白质之间的联系，其中用entrez ID表述基因

2. 表型网络构建：inner_phenotype_similarity.mat文件中的5080\*5080矩阵描述了表型之间两两的相似度（直接用matlab读取）

   **~2~ inner_phenotype_similarity.mat文件中的5080\*5080矩阵对应的OMIM号存储在inter_phenotype_list.txt中,OMIM号对应的疾病名称存储在inter phenotype list explain.txt中*

3. 表型-基因联系：inner_phenotype_gene_relation.txt 左侧第一列数字是phenotype的OMIM号，第二列及右边所有列是该phenotype对应的已知基因序号

   **~3~inner_phenotype_gene_relation.txt中用序号表示基因，而不是entrez ID，需要利用inter_gene_list.txt进行转换，inter_gene_list.txt中存储的数字是entrez ID，元素的索引是改entrez ID对应基因的序号*

​    **~3~inner_phenotype_gene_relation.txt 中共有1126个phenotype，是因为不是所有的 phenotype都    有已知致病基因*

4. inner_ppi.txt 是cipher_main.m程序中的输入，是未扩展前的ppi；

​       ppi_processed_final.txt是扩展后的ppi；

​      cipher_main.m是改进过的cipher算法程序

## 程序解释

```matlab
%% Initialization
clc
clear
%读取phenotype-gene relation中基因序号对应的Entrez ID
gene_list = load('inter_gene_list.txt');

%读取inner_phenotype_similarity.mat中矩阵序号对应的表型OMIM
pheno_list = load('inter_phenotype_list.txt');

inner_pheno_pheno_similarity = load('inner_phenotype_similarity.mat');
fname = fieldnames(inner_pheno_pheno_similarity);
inner_pheno_pheno_similarity = inner_pheno_pheno_similarity.(fname{1});
phenoNum = length(inner_pheno_pheno_similarity);%表型相似度矩阵中的表型数目

clear fname
```

STEP1：先将表格读入matlab中，得到如下数据：

|  gene_list  |  pheno_list   | inner_pheno_pheno_similarity |
| :---------: | :-----------: | :--------------------------: |
| 基因entrez ID | 表型矩阵中对应的OMIM号 |          表型-表型相似度矩阵          |





```matlab
%% Build PPI Network
disp('Computing The Shortest Distance Between Protein Nodes...');
tic

% Load
ppi_matrix_data = load('inner_ppi.txt');%PPI矩阵，用的是entrez ID 表示基因
geneNum = max(max(ppi_matrix_data));%基因数目

%构建稀疏矩阵[sparse函数语法：S = sparse(i,j,v) generates a sparse matrix S from the triplets i, j, and v such that S(i(k),j(k)) = v(k). S = sparse(i,j,v,m,n) specifies the size of S as m-by-n.]此处因为基因间只要有联系就设置值为1，所以v=1；
ppi_matrix = sparse(ppi_matrix_data(:,1),ppi_matrix_data(:,2),1, geneNum, geneNum);

ppi_matrix = ppi_matrix + ppi_matrix';%因为我们要得到的是一个对称矩阵

% Compute PPI_Shortest_Dist
Protein_Shortest_Distance = graphallshortestpaths(ppi_matrix,'Directed','false');

toc
clear ppi_matrix_data
```

STEP2: 构建PPI网络，并计算基因间的距离L~gg'~. 读入ppi.txt得到基因间联系：ppi_matrix_data，因为并不是所有基因两两间都有联系，所以转换为矩阵形式后有大量的零值，因此构建稀疏矩阵以方便存储和调用。然后利用graphallshortestpaths计算；

> `[*dist*] = graphallshortestpaths(*G*)` <u>finds the shortest paths between every pair of nodes in the graph represented by matrix *G*</u>, using Johnson's algorithm. Input *G* is an N-by-N sparse matrix that represents a graph. Nonzero entries in matrix *G* represent the weights of the edges.
>
>   Output *dist* is an N-by-N matrix where `*dist*(S,T)` is the distance of the shortest path from source node `S` to target node `T`. Elements in the diagonal of this matrix are always `0`, indicating the source node and target node are the same. A `0` not in the diagonal indicates that the distance between the source node and target node is `0`. An `Inf` indicates there is no path between the source node and the target node.



```matlab
%% Calculate Phenotype Gene Closeness
disp('Calculate Phenotype Gene Closeness...');
tic
pheno_gene_relation = cell(1,phenoNum);

% Load
fid = fopen('inner_phenotype_gene_relation.txt'); %Get file id——得到文件句柄值
line = fgetl(fid);
while ischar(line) && ~isempty(line)
    chArray = regexp(line,'\t','split'); %Get target array for this drug
    phenoIndex = str2double(chArray(1));
    for k=1:phenoNum
        if(pheno_list(k)==phenoIndex)
            phenoIndex=k;
	    break;
	end
    end
    gNum = size(chArray,2) - 2;
    Array = zeros(1, gNum);
    for i = 1:gNum
        Array(i) = str2double(chArray{i+1});
    end
    pheno_gene_relation{phenoIndex} = Array;
    
    % Next line
    line = fgetl(fid);
end
fclose(fid);
clear fid line chArray Array GeneNum ans i

% Compute Gene2Phenotype Closeness
gene2phenotype_closeness = zeros(geneNum, phenoNum);
for phenoIndex = 1:phenoNum
    Array = pheno_gene_relation{phenoIndex};
    if isempty(Array)
        gene2phenotype_closeness(:, phenoIndex) = 0;
    else
        for geneIndex = 1:geneNum
            gene2phenotype_closeness(geneIndex, phenoIndex) = sum(exp(-(Protein_Shortest_Distance(geneIndex,Array)).^2));
        end
    end
end

toc
clear geneIndex phenoIndex ans array gNum
```

> ***fopen函数***
>
> 1）打开文件
>
> 在读写文件之前，必须先用fopen函数打开或创建文件，并指定对该文件进行的操作方式。fopen函数的调用格式为：
> ```fid=fopen（文件名，‘打开方式’）```
> 说明：其中fid用于**<u>存储文件句柄值</u>**，如果返回的句柄值大于0，则说明文件打开成功。文件名用字符串形式，表示待打开的数据文件。常见的打开方式如下：
>  ‘r’：只读方式打开文件（默认的方式），该文件必须已存在。
>  ‘r+’：读写方式打开文件，打开后先读后写。该文件必须已存在。
>  ‘w’：打开后写入数据。该文件已存在则更新；不存在则创建。
>  ‘w+’：读写方式打开文件。先读后写。该文件已存在则更新；不存在则创建。
>  ‘a’：在打开的文件末端添加数据。文件不存在则创建。
>  ‘a+’：打开文件后，先读入数据再添加数据。文件不存在则创建。
> 另外，在这些字符串后添加一个“t”，如‘rt’或‘wt+’，则将该文件以文本方式打开；如果添加的是“b”，则以二进制格式打开，这也是fopen函数默认的打开方式。
>
> 2）关闭文件
> 文件在进行完读、写等操作后，应及时关闭，以免数据丢失。关闭文件用fclose函数，调用格式为：
> sta＝fclose(fid)
> 说明：该函数关闭fid所表示的文件。sta表示关闭文件操作的返回代码，若关闭成功，返回0，否则返回-1。如果要关闭所有已打开的文件用fclose(‘all’)。
>

> ***fgetl函数***
>
> ```tline=fgetl(fid)  %从文件中读取行，删除文件换行符```
>
> 返回由文件标识符fid指示的文件的下一行。如果fgetl遇到文件结束指示符，则返回-1。对于fid的完整描述请参考fopen函数。fgetl函数常用于含有文件换行符的文件。
>
> ***fgets函数***
>
> tline = fgets(fileID)
>
> 从文件中读取行，保留换行符 (换行符和回车符)
>
> 从文件中读取行，保留*换*行符 读取指定的文件的下一行，包括换行符。 fileid是一个整数文件标识符从*fopen*获得。 tline是一个文本字符串，除非该行只包含结束的文件标记。在这种情况下，tline是数字值-1。与fgets读取字符的编码方案使用与该文件相关联。要指定的编码方案，使用fopen。
>
> 例如：mm.txt文件内容
>
> 1 2 2 3
> 4 5 6
> 2 5 6 8
> 265
> 3
>
> * 利用 fgetl()读入时结果：
>
> \>>c=fgetl(fid)
>
> c =
>
> 1 2 2 3
>
> * 利用fgets()读入结果
>
> \>> a=fgets(fid)
>
> a =
>
> 1 2 2 3
>
> \>> whos c
> Name Size Bytes Class
>
> c 1x7 14 char array
>
> Grand total is 7 elements using 14 bytes
>
> \>> whos a
> Name Size Bytes Class
>
> a 1x9 18 char array
>
> Grand total is 9 elements using 18 bytes
>
> 这就是fgets()为什么比fgetl()多两个字符的原因，在每个换行的时候都会有换行符和回车符。