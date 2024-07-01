clear;
clc;

Data1=readmatrix('UScompHistData18Aug2018-17Aug2023.csv');

Data=Data1(2:927,1:55);

ReturnMatrix=tick2ret(Data);

T=height(ReturnMatrix);
n=width(ReturnMatrix);

%it a objective function which inludes first 749 rows of d_1 ,d_2
%,....d_749 and then 55 of x ,x_1,x_2,......x_55
%and values of d_i(i=1 to 750 is 1/750) and for x_j(j=1 to 55 is 0)
f1=1/T*(ones(T,1));
f2=zeros(n,1);
f=cat(1,f1,f2)';

%AvgReturn variable conatins an average of individual company
AvgReturn=(1/T)*sum(ReturnMatrix,1);

%A which contains fiirst 749 +ve , then 749 -ve and 1 for meu not
Positive_A=zeros(height(ReturnMatrix) , height(ReturnMatrix)+width(ReturnMatrix));
Negative_A=Positive_A;
Positive_A(1 : T, 1 : T) = (-1) * eye(T);
Negative_A(1 : T, 1 : T) = (-1) * eye(T);
for i = 1 : T
    for j = T+1 : T+n
        Positive_A(i,j)=ReturnMatrix(i,j-T)-AvgReturn(j-T);
        Negative_A(i,j)=(-1)*ReturnMatrix(i,j-T)+AvgReturn(j-T);
    end
end    
LastCon=cat(2, zeros(1, T), (-1) * AvgReturn);
A=cat(1,Positive_A,Negative_A,LastCon);
b=zeros(1,2*T+1)';
b(height(b),1)=-0.001;
Aeq=cat(2 , zeros(1,T) , ones(1,n));
beq=1;
lb=zeros(1,T+n);
ub=[];
[X,Z]=linprog(f,A,b,Aeq,beq,lb,ub);
X_values=X(T+1:T+n);
Z;