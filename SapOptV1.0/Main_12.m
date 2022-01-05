% *****************************************************************************
% 
%  Main_12.m
% 
%  Description:
%      This program is used to optimize the disk in aero engine.
% 
% 
%  All the codes had been written by Zhizhen Dong in 2020-2021
% 
% *****************************************************************************

clc;
clear;
close all;
%26.229s

%% 导入数据、初始化
%周向网格层数,单元层数（输入！！！）
zn=5;

%网格的中间节点数（输入！！！）
mn=3; 

% 导入数据///可改为根据节点编号寻找
load couplednodes.txt         %导入耦合节点编号（FROM APDL）
left=[863;789;couplednodes(:,1);907;902];       %需要观察补充
right=[774;777;couplednodes(:,mn+2);882;876];  %需要观察补充

% 节点关联
left_coup_index=[3 length(left)-2]';         %起始节点、终止节点索引
right_coup_index=[3 length(right)-2]';        %起始节点、终止节点索引

%最大限制应力
limit_body_eqv=270;
limit_root_eqv=300;

%% 导入全部节点信息（编号、坐标）换算
nodeinitial=xlsread('nodedata.xlsx','data3','A1:D1016'); %导入左边界节点         、坐标
nn=nodeinitial(:,1);
x=nodeinitial(:,2);
y=nodeinitial(:,3);
z=nodeinitial(:,4);

[th,r,z] = cart2pol(x,z,y);

nodetrans=[r,z,nn,th];  %节点、单元全部信息整合

for i=1:1:length(left(:,1))              %已知节点编号查left网格索引
    [left_index(i,:),]=find(nodetrans(:,3)==left(i,1));
end

for i=1:1:length(right(:,1))             %已知节点编号查right网格索引
    [right_index(i,:),]=find(nodetrans(:,3)==right(i,1));
end

%% 关联2D平面待更新耦合节点与3D体待更新耦合节点
temp=1;
coupled_2D_index=zeros(size(couplednodes,1)*size(couplednodes,2),1); %初始化全部耦合节点在APDL网格中的索引
for i=1:size(couplednodes,2)
    for j=1:size(couplednodes,1)
        coupled_2D_index(temp)=find(nodetrans(:,3)==couplednodes(j,i));
        temp=temp+1;
    end
end

coupled_2D_rz=[nodetrans(coupled_2D_index,1),nodetrans(coupled_2D_index,2)];  

coupled_3D_index=zeros(size(couplednodes,1)*size(couplednodes,2),zn+1);   %初始化
coupled_3D=zeros(size(couplednodes,1)*size(couplednodes,2),zn+1);   %初始化

for i=1:size(couplednodes,1)*size(couplednodes,2)
    tt=1;
    temp_r_index=find(nodetrans(:,1)<(coupled_2D_rz(i,1)+0.005)&(coupled_2D_rz(i,1)-0.005)<nodetrans(:,1));
    temp_z_index=find(nodetrans(:,2)<(coupled_2D_rz(i,2)+0.001)&(coupled_2D_rz(i,2)-0.001)<nodetrans(:,2));
    for j=1:length(temp_r_index)   %通用性提高
    	if find(temp_r_index(j)==temp_z_index(:)) %    && find(temp_rz_index(j)==temp_z_index(:))
            temp_index(tt,1)=temp_r_index(j);
            tt=tt+1;
        end
    end

    for j=1:length(temp_index)
        coupled_3D_index(i,:)=temp_index;      %耦合节点索引
    end
end

for i=1:size(coupled_3D_index,1)
    for j=1:size(coupled_3D_index,2)
        coupled_3D(i,j)=nodetrans(coupled_3D_index(i,j),3);         %耦合节点号
    end
end

%还有问题 排序问题！！！
coupled_3D=[coupled_3D(:,6),coupled_3D(:,5),coupled_3D(:,1),coupled_3D(:,2),coupled_3D(:,3),coupled_3D(:,4)];
coupled_3D_index=[coupled_3D_index(:,6),coupled_3D_index(:,5),coupled_3D_index(:,1),coupled_3D_index(:,2),coupled_3D_index(:,3),coupled_3D_index(:,4)];

%% 筛选除待更新节点外的固定节点
nodefixed=nodeinitial;  
nodefixed(coupled_3D_index,:)=[]; 

%% 初始化
left_min_z=nodetrans(couplednodes(:,1),2)-[0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5]-[0.5;1;1;1;1;1;1;1;1;1;1;1;1;0.5];      %参数搜索范围
left_max_z=nodetrans(couplednodes(:,1),2)+[2;3;4;4.8;4.8;4.8;4.8;4.8;4.8;4.8;4.8;4;3;2]+[0.5;1;1.2;1.2;1.5;3;3;2.5;2;1.5;1.2;1.2;1;0.5];

right_min_z=nodetrans(couplednodes(:,mn+2),2)-[2;3;4;4.8;4.8;4.8;4.8;4.8;4.8;4.8;4.8;4;3;2]-[0.5;1;1.2;1.2;1.5;2.5;2.5;3;3;1.5;1.2;1.2;1;0.5];      %参数搜索范围
right_max_z=nodetrans(couplednodes(:,mn+2),2)+[0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5]+[0.5;1;1;1;1;1;1;1;1;1;1;1;1;0.5];

MinX=[left_min_z;right_min_z];
MaxX=[left_max_z;right_max_z];

Vmax=1;              %限定速度的范围
Vmin=-1;

Size=50;      
CodeL=length(left_min_z)+length(right_min_z);         %参数个数

c1=1.9;c2=2.0;c3=0.1;c4=0.15;   %学习因子
wmax=0.9;wmin=0.4; 
G=500;         %最大迭代次数

for i=1:G
    w(i)=wmax-((wmax-wmin)/G)*i;  
end

for i=1:1:CodeL/2
    X(:,i)=nodetrans(couplednodes(i,1),2)+(left_max_z(i)-left_min_z(i))*rand(Size,1);              
    X(:,CodeL/2+i)=nodetrans(couplednodes(i,mn+2),2)-(right_max_z(i)-right_min_z(i))*rand(Size,1);                                    
    v(:,i)=Vmin+(Vmax-Vmin)*rand(Size,1);                                          
    v(:,CodeL/2+i)=Vmin+(Vmax-Vmin)*rand(Size,1);
end

%%
for i=1:1:Size
	ssaleft=ssa([nodetrans(863,2),nodetrans(789,2),X(i,1:CodeL/2),nodetrans(907,2),nodetrans(902,2)]); %光顺
	X(i,1:CodeL/2)=ssaleft(3:length(ssaleft)-2);
	ssaright=ssa([nodetrans(774,2),nodetrans(777,2),X(i,CodeL/2+1:CodeL),nodetrans(882,2),nodetrans(876,2)]);
	X(i,CodeL/2+1:CodeL)=ssaright(3:length(ssaright)-2);
end

%% 节点映射
mapping

%% 节点位置输出
export

%% ANSYS求解
system('"D:\Tools\ANSYS Inc\v202\ansys\bin\winx64\ANSYS202" -b -p ane3fl -i C:\Users\DongZZ\Desktop\optResearch_1st\mixedcode\Main_07.f -o D:\DZZ.out') 

%% 最优解集初始化
load matapdlrst.txt         %导入求解结果（FROM APDL）
body_eqv=matapdlrst(:,1);
root_eqv=matapdlrst(:,2);
disc_volu=matapdlrst(:,3);
score=disc_volu*0.06+0.5*(body_eqv-limit_body_eqv).^2+0.5*(root_eqv-limit_root_eqv).^2;

%  个体最佳初始化
Xhb=X;             %初始化当前每组粒子最佳位置Xl(i,:)为初始位置X(i,:)  historybest
scorehb=score;

%  群体最佳初始化
system('"D:\Tools\ANSYS Inc\v202\ansys\bin\winx64\ANSYS202" -b -p ane3fl -i C:\Users\DongZZ\Desktop\optResearch_1st\mixedcode\precal.f -o D:\DZZ.out')
load matapdlini.txt

BestS=[nodetrans(couplednodes(:,1),2);nodetrans(couplednodes(:,mn+2),2)]';                %初始化初始边界及对应计算结果为全局最优
scorehgb=matapdlini(3)*0.06+0.5*(matapdlini(1)-limit_body_eqv)^2+0.5*(matapdlini(2)-limit_root_eqv)^2;


%% 粒子群更新主程序
BestJ=zeros(G,1);
optPath=zeros(G,CodeL);

disPN=scorehgb-1;    
hisBestScore=scorehgb;
hwait=waitbar(0,'总进度>>>>>>>>'); 
for kg=1:1:G             
    BestJ(kg)=vpa(scorehgb/hisBestScore,4);
    optPath(kg,:)=BestS;
    for i=1:Size         
        v(i,:)=w(kg)*v(i,:)+c1*rand*(Xhb(i,:)-X(i,:))+c2*rand*(BestS-X(i,:));%+c3*rand*([left_max_z;right_min_z]'-X(i,:));%c3*rand*([(X(i,1:CodeL/2)+X(i,CodeL/2+1:CodeL))./2,(X(i,1:CodeL/2)+X(i,CodeL/2+1:CodeL))./2]-X(i,:));%+c4*rand*([nodetrans(couplednodes(:,1),2);nodetrans(couplednodes(:,mn+2),2)]'-X(i,:));       %更新当前组粒子速度
        v(i,:)=v(i,:)+[0.01*ones(1,CodeL/2),-0.01*ones(1,CodeL/2)]*rand;
        for j=1:CodeL               %规范每组各参数粒子速度
            if v(i,j)<Vmin
                v(i,j)=Vmin;
            elseif v(i,j)>Vmax
                v(i,j)=Vmax;
            end
        end
        X(i,:)=X(i,:)+v(i,:);       %更新当前组粒子位置
        for j=1:CodeL               %规范每组各参数粒子位置
            if X(i,j)<MinX(j)
                X(i,j)=MinX(j);
            elseif X(i,j)>MaxX(j)
                X(i,j)=MaxX(j);
            end
        end 

        %随机数跳出局部最优解
        if kg>50 && kg<100 && rand>0.95              %5%的概率触发     
            k=ceil(CodeL*rand);      %k的值为1-4(4个参数里面随机变异)  
            if k<CodeL/2+1
                X(i,k)=left_min_z(k)+(left_max_z(k)-left_min_z(k))*rand;       %当前粒子组，第k个参数变异为0-2之间的位置 
            else
                X(i,k)=right_min_z(k-CodeL/2)+(right_max_z(k-CodeL/2)-right_min_z(k-CodeL/2))*rand;
            end
        elseif kg>100 && kg<150 && rand>0.9              %10%的概率触发            
            k=ceil(CodeL*rand);      %k的值为1-4(4个参数里面随机变异)  
            if k<CodeL/2+1
                X(i,k)=left_min_z(k)+(left_max_z(k)-left_min_z(k))*rand;       %当前粒子组，第k个参数变异为0-2之间的位置 
            else
                X(i,k)=right_min_z(k-CodeL/2)+(right_max_z(k-CodeL/2)-right_min_z(k-CodeL/2))*rand;
            end
        elseif kg>150 && rand>0.85              %15%的概率触发            
            k=ceil(CodeL*rand);      %k的值为1-4(4个参数里面随机变异)  
            if k<CodeL/2+1
                X(i,k)=left_min_z(k)+(left_max_z(k)-left_min_z(k))*rand;       %当前粒子组，第k个参数变异为0-2之间的位置 
            else
                X(i,k)=right_min_z(k-CodeL/2)+(right_max_z(k-CodeL/2)-right_min_z(k-CodeL/2))*rand;
            end         
        end
    end   
        
	if rem(kg,15)==1 && kg<G/2
        for i=1:1:Size
            ssaleft=ssa([nodetrans(863,2),nodetrans(789,2),X(i,1:CodeL/2),nodetrans(907,2),nodetrans(902,2)]); 
            X(i,1:CodeL/2)=ssaleft(3:length(ssaleft)-2);
            ssaright=ssa([nodetrans(774,2),nodetrans(777,2),X(i,CodeL/2+1:CodeL),nodetrans(882,2),nodetrans(876,2)]);
            X(i,CodeL/2+1:CodeL)=ssaright(3:length(ssaright)-2);
        end
	end   
                
    %输出粒子位置
    mapping
    export

    %ANSYS求解
    system('"D:\Tools\ANSYS Inc\v202\ansys\bin\winx64\ANSYS202" -b -p ane3fl -i C:\Users\DongZZ\Desktop\optResearch_1st\mixedcode\Main_07.f -o D:\DZZ.out')
    
	%结果导入
	load matapdlrst.txt         %导入求解结果（FROM APDL）
	body_eqv=matapdlrst(:,1);
	root_eqv=matapdlrst(:,2);
	disc_volu=matapdlrst(:,3);
	score=disc_volu*0.06+0.5*(body_eqv-limit_body_eqv).^2+0.5*(root_eqv-limit_root_eqv).^2;

    for i=1:Size         %粒子组数
        if scorehb(i)>score(i)  
            Xhb(i,:)=X(i,:);
            scorehb(i)=score(i);          %策略加权
            if scorehgb>score(i)                                      
                BestS=X(i,:);  
                scorehgb=score(i);      %策略加权
            end
        end
    end  
      
    %进度条更新
    str=['主程序运行中',num2str(kg/G*100),'%']; %用于显示  
	waitbar(kg/G,hwait,str);  
	%pause(0.05); %pause(n)函数是程序停止n秒后继续
 
end
BestJ
optPath
save opt.mat