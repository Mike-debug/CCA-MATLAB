%����ʼ����������ڴ沢����
clear;clc;

%��λ���ݴ�ŵ��ļ�
XFile='XFile.txt';
YFile='YFile.txt';
%�������ݵ�Ҫ��ά�ȴ���1�������ݼ�X��Y������ͬ����������
%��������
X=importdata(XFile);
Y=importdata(YFile);

%{
X,Y�����ڶ���֮��������������һ������������������ͬһά���ϸ���������ֵ��
�ں����ļ����л��X��Yת�ã���������matlab�⺯��zscore����Ҫ��������δת�õľ�����б�׼���ļ���
%}
X=zscore(X);                                                               %��׼��X����������ֵΪ0������Ϊ1
Y=zscore(Y);                                                               %��׼��Y����������ֵΪ0������Ϊ1

%ת��X��Y
X_=X';
Y_=Y';
%{
%�����׼�������ݾ���Ĵ�С
Size=size(X_);
%��ʾ��׼�������ݵ�
for i=1:Size(1,1)
    U1=X_(i,:);
    U1=U1';
    V1=Y_(i,:);
    V1=V1';
    figure(i);
    plot(V1,U1,'*');
end
%}

Sxx=X_*X_';%X_�ķ���
Syy=Y_*Y_';%Y_�ķ���
Sxy=X_*Y_';%X_��Y_��Э����
gama=0.0;%����ϵ����
manzhix=length(Sxx);%Sxx������
manzhiy=length(Syy);%Syy������
%{
Sxx��Syy����һ�������ȣ�������������һ����Ϊ�趨������ϵ��
��������Ϊ����ϵ����Ӧ������С��������1e-7Ϊ��С��λ��ʼ����
%}
while(rank(Sxx)~=manzhix || rank(Syy)~=manzhiy )
    Sxx=Sxx+0.000001*eye(manzhix);
    Syy=Syy+0.000001*eye(manzhiy);
    gama=gama+1;
end

M = (Sxx^ - 0.5)*Sxy*(Syy^ - 0.5);
[U, S, V]=svd(M);%����ֵ�ֽ�

%ѭ����ʾ��������ֵ�µ����Ա任�Լ��ɸ����Ա任����õ������ϵ����
for i=1:min(size(S))
    u = U(:,i);                                                                 %�������ֵ������������
    v = V(:,i);                                                                 %�������ֵ������������
    a=(Sxx^-0.5)*u;                                                             %X_ͶӰ����
    b=(Syy^-0.5)*v;                                                             %Y_ͶӰ����
    X1=a'*X_;
    Y1=b'*Y_;

    if var(X1)<=1e-10||var(Y1)<=1e-10
        rou=1;
    else
        rou=mean(X1.*Y1)-mean(X1)*mean(Y1);                                                          %������ϵ��
        rou=rou/std(X1,1,2);
        rou=rou/std(Y1,1,2);
    end
%{
&�˶δ���ԭ�������׼��������ݾ���ͱ�׼�������Ա任������ݾ���
    disp('X_=');
    disp(X_);
    disp('Y_=');
    disp(Y_);
    disp('X1=');
    disp(X1);
    disp('Y1=');
    disp(Y1);
%}
%��ʾ�ڼ������Ա任
    fprintf('=====================================��%d�����Ա任=====================================\n',i);
%����������Ա任
    disp('���Ա任Ϊ��');
    fprintf('a = \n');
    disp(a);
    fprintf('b = \n');
    disp(b);
%������ϵ��
    fprintf('rou = \n');
    disp(rou);
%�������ϵ��
    fprintf('gama = \n');
    disp(gama);
%�������
    fprintf('\n\n\n');
    
    
%��ͼ
    X2=a'*X_;%�������Ա任��ı�׼������
    Y2=b'*Y_;
    figure(i)%ͼi;
    plot(X2,Y2,'*');%��ɢ����ͼ
end