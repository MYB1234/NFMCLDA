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