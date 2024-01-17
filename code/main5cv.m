clc
clear
%%%load data
interactionA=xlsread('...\data\05-lncRNA-disease');
x=gKernel(240,interactionA);
y=gKernel(412,interactionA');
sl1=xlsread('...\data\04-lncRNA-lncRNA');
sl=(x+sl1)*0.5;
L=transition_matrix_network_propagation1(sl);

sd1=xlsread('...\data\07-disease-disease');
sd=(y+sd1)*0.5;
D=transition_matrix_network_propagation1(sd);

interactionC=xlsread('...\data\06-miRNA-disease');
sm=gKernel(495,interactionC);
M=transition_matrix_network_propagation1(sm);

interactionB=xlsread('...\data\08-lncRNA-miRNA');

A_ori=interactionA;

score=thrrw(L,D,interactionA,interactionB,interactionC,M, 1,3, 1,0.3,0.1,0.9);
T=score;
trIndex = double(T ~= 0);
[WW,iter] = MC(0.9,0.6, 100, T, trIndex, 0.001, 300, 0,1);

seed = 1;
%%%%
CV=5;    %%We take 5-fold crossvalidation 
rand('state',seed);
[row,col]=size(A_ori);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[row_index,col_index]=find(A_ori==1); %%find all the elments that its value is 1
link_num=sum(A_ori(:));  %% caculate the numbers of the interaction
rand('state',seed);
random_index=randperm(link_num);
size_of_CV=round(link_num/CV);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result=zeros(1,7);
for k=1:CV  %%%
    fprintf('begin to implement the cross validation:round =%d/%d\n', k, CV);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (k~=CV) %% We allocate all the interaction elements into 5 parts
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
    else
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):end));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):end));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    train_set=A_ori;
    test_link_num=size(test_row_index,1);
    for i=1:test_link_num %% let interaction elements' value equal to 0 in the test matrix
        train_set(test_row_index(i),test_col_index(i))=0;  %
    end
    
    score=thrrw(L,D,train_set,interactionB,interactionC,M, 1, 3, 1,0.3,0.1,0.9);
    T=score;
    trIndex = double(T ~= 0);
    [WW,iter] = MC(0.9,0.6, 100, T, trIndex, 0.001, 300, 0,1);
    result=result+model_evaluate(A_ori,WW,train_set);
    result/k
end
result=result/CV;
disp('aupr,auc,sen,spec,precision,accuracy,f1')
result
