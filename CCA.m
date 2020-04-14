%程序开始，清除所有内存并清屏
clear;clc;

%定位数据存放的文件
XFile='XFile.txt';
YFile='YFile.txt';
%对于数据的要求，维度大于1，且数据集X和Y具有相同的数据组数
%导入数据
X=importdata(XFile);
Y=importdata(YFile);

%{
X,Y矩阵在读入之初，行向量代表一个样本，列向量代表同一维度上各个样本的值。
在后续的计算中会把X和Y转置，出于利用matlab库函数zscore的需要，这里用未转置的矩阵进行标准化的计算
%}
X=zscore(X);                                                               %标准化X，列向量均值为0，方差为1
Y=zscore(Y);                                                               %标准化Y，列向量均值为0，方差为1

%转置X和Y
X_=X';
Y_=Y';
%{
%计算标准化后数据矩阵的大小
Size=size(X_);
%显示标准化后数据的
for i=1:Size(1,1)
    U1=X_(i,:);
    U1=U1';
    V1=Y_(i,:);
    V1=V1';
    figure(i);
    plot(V1,U1,'*');
end
%}

Sxx=X_*X_';%X_的方差
Syy=Y_*Y_';%Y_的方差
Sxy=X_*Y_';%X_和Y_的协方差
gama=0.0;%正则化系数γ
manzhix=length(Sxx);%Sxx的满秩
manzhiy=length(Syy);%Syy的满秩
%{
Sxx和Syy中有一个不满秩，则两个都加上一个人为设定的正则化系数
本程序认为正则化系数γ应当尽量小，所以以1e-7为最小单位开始加起
%}
while(rank(Sxx)~=manzhix || rank(Syy)~=manzhiy )
    Sxx=Sxx+0.000001*eye(manzhix);
    Syy=Syy+0.000001*eye(manzhiy);
    gama=gama+1;
end

M = (Sxx^ - 0.5)*Sxy*(Syy^ - 0.5);
[U, S, V]=svd(M);%奇异值分解

%循环显示各个奇异值下的线性变换以及由各线性变换计算得到的相关系数ρ
for i=1:min(size(S))
    u = U(:,i);                                                                 %最大奇异值的左奇异向量
    v = V(:,i);                                                                 %最大奇异值的右奇异向量
    a=(Sxx^-0.5)*u;                                                             %X_投影向量
    b=(Syy^-0.5)*v;                                                             %Y_投影向量
    X1=a'*X_;
    Y1=b'*Y_;

    if var(X1)<=1e-10||var(Y1)<=1e-10
        rou=1;
    else
        rou=mean(X1.*Y1)-mean(X1)*mean(Y1);                                                          %最大相关系数
        rou=rou/std(X1,1,2);
        rou=rou/std(Y1,1,2);
    end
%{
&此段代码原本输出标准化后的数据矩阵和标准化和线性变换后的数据矩阵
    disp('X_=');
    disp(X_);
    disp('Y_=');
    disp(Y_);
    disp('X1=');
    disp(X1);
    disp('Y1=');
    disp(Y1);
%}
%提示第几个线性变换
    fprintf('=====================================第%d个线性变换=====================================\n',i);
%输出两个线性变换
    disp('线性变换为：');
    fprintf('a = \n');
    disp(a);
    fprintf('b = \n');
    disp(b);
%输出相关系数
    fprintf('rou = \n');
    disp(rou);
%输出正则化系数
    fprintf('gama = \n');
    disp(gama);
%输出空行
    fprintf('\n\n\n');
    
    
%绘图
    X2=a'*X_;%经过线性变换后的标准化数据
    Y2=b'*Y_;
    figure(i)%图i;
    plot(X2,Y2,'*');%离散点阵图
end