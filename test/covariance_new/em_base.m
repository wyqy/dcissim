%% https://zhuanlan.zhihu.com/p/438941218
%% 基本使用
PntSet1=mvnrnd([2 3],[1 0;0 2],500);
PntSet2=mvnrnd([6 7],[1 0;0 2],500);
PntSet3=mvnrnd([14 9],[1 0;0 1],500);
PntSet=[PntSet1;PntSet2;PntSet3];

% 构造GMM模型
tic
[Mu,Sigma,Pi,Class]=gaussKMeans(PntSet,3,'dis');
toc

%% 进阶使用
% 构造概率密度函数
func=getGaussFunc(Mu,Sigma,Pi);

% 绘制概率密度图像
figure('Units','normalized','Position',[.3,.2,.6,.65])
[X1,X2]=meshgrid(0:.4:16,0:.4:12);
surf(X1,X2,func(X1,X2),'LineWidth',1)
%修饰一下
ax=gca;hold(ax,'on');
ax.XLim=[0,16];
ax.YLim=[0,12];
ax.LineWidth=2;
ax.Box='on';
ax.TickDir='in';
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.ZMinorTick='on';
ax.XColor=[.3,.3,.3];
ax.YColor=[.3,.3,.3];
ax.ZColor=[.3,.3,.3];
ax.FontWeight='bold';
ax.FontName='Cambria';
ax.FontSize=13;
ax.GridLineStyle='--';

%% 绘图
% 绘制散点图
figure('Units','normalized','Position',[.3,.2,.5,.65])
ax=gca;hold(ax,'on');
colorList=[0.4  0.76 0.65
           0.99 0.55 0.38 
           0.55 0.63 0.80];
for i=1:3
    scatter(PntSet(Class==i,1),PntSet(Class==i,2),180,'filled',...
        'LineWidth',2.2,'MarkerEdgeColor',[1 1 1]*.3,'MarkerFaceColor',colorList(i,:));
end
% 绘制置信椭圆       
for i=1:3
    [X,Y]=getEllipse(Mu(i,:),Sigma(:,:,i),9.21,100);
    fill(X,Y,colorList(i,:),'EdgeColor',colorList(i,:).*.5,...
        'LineWidth',3,'FaceAlpha',.2)
    
end     
legend('pointSet1','pointSet2','pointSet3','conf.ellipse1','conf.ellipse2','conf.ellipse3')
%修饰一下
ax.XLim=[-2,18];
ax.YLim=[-2,14];
ax.LineWidth=2;
ax.Box='on';
ax.TickDir='in';
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.XGrid='on';
ax.YGrid='on';
ax.GridLineStyle='--';
ax.XColor=[.3,.3,.3];
ax.YColor=[.3,.3,.3];
ax.FontWeight='bold';
ax.FontName='Cambria';
ax.FontSize=15;
%ax.Color=[0.9,0.9,0.9];



%% Gauss-EM函数
function [Mu,Sigma,Pi,Class]=gaussKMeans(pntSet,K,initM)
% @author:slandarer
% ===============================================================
% pntSet  | NxD数组   | 点坐标集                                |
% K       | 数值      | 划分堆数量                              |
% --------+-----------+-----------------------------------------+
% Mu      | KxD数组   | 每一行为一类的坐标中心                  |
% Sigma   | DxDxK数组 | 每一层为一类的协方差矩阵                |
% Pi      | Kx1列向量 | 每一个数值为一类的权重(占比)            |
% Class   | Nx1列向量 | 每一个数值为每一个元素的标签(属于哪一类)|
% --------+-----------+-----------------------------------------+


[N,D]=size(pntSet); % N:元素个数 | D:维数

% 初始化数据===============================================================
if nargin<3
    initM='random';
end
switch initM
    case 'random' % 随机取初始值
        [~,tIndex]=sort(rand(N,1));tIndex=tIndex(1:K);
        Mu=pntSet(tIndex,:);

    case 'dis'    % 依据各维度的最大最小值构建方向向量
                  % 并依据该方向向量均匀取点作为初始中心       
        tMin=min(pntSet);
        tMax=max(pntSet);
        Mu=linspace(0,1,K)'*(tMax-tMin)+repmat(tMin,K,1);

    % case '依据个人需求自行添加'  
    % ... ...
    % ... ...     
end

% 一开始设置每一类有相同协方差矩阵和权重
Sigma(:,:,1:K)=repmat(cov(pntSet),[1,1,K]);
Pi(1:K,1)=(1/K);

% latest coefficient:上一轮的参数
LMu=Mu;        
LPi=Pi;
LSigma=Sigma;

turn=0; %轮次

% GMM/gauss_k_means主要部分================================================
while true
    
    % 计算所有点作为第k类成员时概率及概率和(不加权重)
    % 此处用了多次转置避免构建NxN大小中间变量矩阵
    % 而将过程中构建的最大矩阵缩小至NxD，显著减少内存消耗
    Psi=zeros(N,K);
    for k=1:K
        Y=pntSet-repmat(Mu(k,:),N,1);
        Psi(:,k)=((2*pi)^(-D/2))*(det(Sigma(:,:,k))^(-1/2))*...
                      exp(-1/2*sum((Y/Sigma(:,:,k)).*Y,2))';    
    end
    
    % 加入权重计算各点属于各类后验概率
    Gamma=Psi.*Pi'./sum(Psi.*Pi',2);
    
    % 大量使用矩阵运算代替循环，提高运行效率
    Mu=Gamma'*pntSet./sum(Gamma,1)';
    for k=1:K
        Y=pntSet-repmat(Mu(k,:),N,1);
        Sigma(:,:,k)=(Y'*(Gamma(:,k).*Y))./sum(Gamma(:,k));
    end
    Pi=(sum(Gamma)/N)';
    [~,Class]=max(Gamma,[],2);

    % 计算均方根误差
    R_Mu=sum((LMu-Mu).^2,'all');
    R_Sigma=sum((LSigma-Sigma).^2,'all');
    R_Pi=sum((LPi-Pi).^2,'all');
    R=sqrt((R_Mu+R_Sigma+R_Pi)/(K*D+D*D*K+K));
    
    % 每隔10轮输出当前收敛情况
    turn=turn+1;
    if mod(turn,10)==0
        disp(' ')
        disp('==================================')
        disp(['第',num2str(turn),'次EM算法参数估计完成'])
        disp('----------------------------------')
        disp(['均方根误差:',num2str(R)])
        disp('当前各类中心点：')
        disp(Mu)
    end
    
    % 循环跳出
    if (R<1e-6)||isnan(R)
        disp(['第',num2str(turn),'次EM算法参数估计完成'])
        if turn>=1e4||isnan(R)
            disp('GMM模型不收敛')
        else
            disp(['GMM模型收敛，参数均方根误差为',num2str(R)])
        end
        break;
    end
    
    LMu=Mu;
    LSigma=Sigma;
    LPi=Pi;
end
end


%% 高维混合分布函数生成函数
function func=getGaussFunc(Mu,Sigma,Pi)
[K,D]=size(Mu);

X{D}=[];
for d=1:D
    X{d}=['x',num2str(d)];
end
X=sym(X);

func=0;
for k=1:K
    tMu=Mu(k,:);
    tSigma=Sigma(:,:,k);   
    tPi=Pi(k);
    tX=X-tMu;   
    func=func+tPi*(1/(2*pi)^(D/2))*(1/det(tSigma)^(1/2))*exp((-1/2)*(tX/tSigma*tX.'));
end

func=matlabFunction(func);
end


%% 椭圆坐标生成函数
function [X,Y]=getEllipse(Mu,Sigma,S,pntNum)
% 置信区间 | 95%:5.991  99%:9.21  90%:4.605
% (X-Mu)*inv(Sigma)*(X-Mu)=S

invSig=inv(Sigma);

[V,D]=eig(invSig);
aa=sqrt(S/D(1));
bb=sqrt(S/D(4));

t=linspace(0,2*pi,pntNum);
XY=V*[aa*cos(t);bb*sin(t)];
X=(XY(1,:)+Mu(1))';
Y=(XY(2,:)+Mu(2))';
end

