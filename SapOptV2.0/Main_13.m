clc;
clear;
close all;

%% 导入数据、初始化
%网格的中间节点数（输入！！！）
mn=4; 

% 导入数据///可改为根据节点编号寻找
load couplednodes.txt         %导入耦合节点编号（FROM APDL）
left=[2956;2947;2948;couplednodes(:,1);2611;2612;2623;2624];       %需要观察补充
right=[2730;2734;2739;couplednodes(:,mn+2);2669;2672;2681;2684];  %需要观察补充

%计算chebyshev点
a=32.5;  %范围
b1=95;
b2=98.8;
n=8;             %阶数-1

left_cheby_r=chebyshev(a,b1,n)';  %R
right_cheby_r=chebyshev(a,b2,n)'; %R

%最大限制应力
limit_body_eqv=375;

%% 导入全部节点信息（编号、坐标）极坐标换算
nodeinitial=xlsread('nodedata.xlsx','datapaper','A1:D3188'); %导入节点、坐标
nn=nodeinitial(:,1);
x=nodeinitial(:,2);
y=nodeinitial(:,3);
z=nodeinitial(:,4);

[th,r,z] = cart2pol(x,z,y);

nodetrans=[r,z,nn,th]; 

for i=1:1:length(left(:,1))              %已知节点编号查left网格索引
    [left_index(i,:),]=find(nodetrans(:,3)==left(i,1));
end

for i=1:1:length(right(:,1))             %已知节点编号查right网格索引
    [right_index(i,:),]=find(nodetrans(:,3)==right(i,1));
end

%%   周向各层角坐标提取
sita_1=ones(length(nodetrans(:,4)),1);
for i=1:length(nodetrans(:,4))
    sita_1(i,1)=round(nodetrans(i,4)*100000)/100000;  
end
sita_2=unique(sita_1);

%倒叙排列，接近实际
sita_3=sita_2(end:-1:1);

%极坐标按角坐标分类
for i=1:length(sita_2)
    tt=1;
    for j=1:length(nodetrans)
        if (nodetrans(j,4)<=sita_2(i)+0.001)&&(nodetrans(j,4)>=sita_2(i)-0.001)
            sita_part(tt,:,i)=nodetrans(j,:);
            tt=tt+1;
        end
    end
end

%倒叙排列，接近实际
sita_part2=sita_part(:,:,end:-1:1);

%周向网格层数,单元层数
zn=length(sita_2)-1;

%%  2D端面内节点关联
%2D端面内节点期望坐标计算
facenodes=xlsread('nodedata.xlsx','dataface','A1:A488'); %导入左边界节点         、坐标
copnds1=xlsread('nodedata.xlsx','core','A1:B48'); 
copnds2=xlsread('nodedata.xlsx','expanded','A1:B6'); 
copnds=[copnds1;copnds2];

deltadis_r=zeros(length(copnds),1);
deltadis_z=zeros(length(copnds),1);
for i=1:length(copnds)   
    discopnds_r=nodetrans(copnds(i,2),1)-nodetrans(copnds(i,1),1);  
    discopnds_z=nodetrans(copnds(i,2),2)-nodetrans(copnds(i,1),2);
    deltadis_r(i)=discopnds_r/(mn+1);
    deltadis_z(i)=discopnds_z/(mn+1);
end

coupled_2D_temp=zeros(length(copnds),mn+2); %初始化全部端面耦合节点在APDL网格中的索引
coupled_2D_ndr=zeros(length(copnds),mn+2); %初始化全部端面耦合节点在APDL网格中的r坐标
coupled_2D_ndz=zeros(length(copnds),mn+2); %初始化全部端面耦合节点在APDL网格中的z坐标

for i=1:length(copnds)
    for j=2:mn+1
        coupled_2D_ndr(i,j)=nodetrans(copnds(i,1),1)+deltadis_r(i)*(j-1);
        coupled_2D_ndz(i,j)=nodetrans(copnds(i,1),2)+deltadis_z(i)*(j-1);      
    end
end
coupled_2D_ndr(:,1)=nodetrans(copnds(:,1),1);
coupled_2D_ndr(:,mn+2)=nodetrans(copnds(:,2),1);
coupled_2D_ndz(:,1)=nodetrans(copnds(:,1),2);
coupled_2D_ndz(:,mn+2)=nodetrans(copnds(:,2),2);

%2D端面内节点索引匹配
deltadis_t=sqrt(deltadis_r.^2+deltadis_z.^2);
trshd=min(deltadis_t)/2;         %比较关键

for i=1:length(copnds)
    for j=1:mn+2
        temp_r_index=find(nodetrans(facenodes(:),1)<(coupled_2D_ndr(i,j)+trshd)&(coupled_2D_ndr(i,j)-trshd)<nodetrans(facenodes(:),1));
        temp_z_index=find(nodetrans(facenodes(:),2)<(coupled_2D_ndz(i,j)+trshd)&(coupled_2D_ndz(i,j)-trshd)<nodetrans(facenodes(:),2));

        for k=1:length(temp_r_index)   %通用性提高
            if find(temp_r_index(k)==temp_z_index(:)) %    && find(temp_rz_index(j)==temp_z_index(:))
                coupled_2D_temp(i,j)=facenodes(temp_r_index(k));
            end
        end
    end
end

tt=1;
coupled_2D_index=zeros(size(coupled_2D_temp,1)*size(coupled_2D_temp,2),1);
for j=1:size(coupled_2D_temp,2)           %数据整理
    for i=1:size(coupled_2D_temp,1)   
        coupled_2D_index(tt)=coupled_2D_temp(i,j);
        tt=tt+1;
    end
end

% %重复索引验证
% kk=mode(coupled_2D_index');
% find(coupled_2D_index==kk)

Size=50;        %粒子组数
coupled_2D_ntr=zeros(length(coupled_2D_index),1);%后续使用
coupled_2D_ntz=zeros(length(coupled_2D_index),Size);%后续使用
coupled_2D_ntr=nodetrans(coupled_2D_index,1);
% for i=1:Size
%     coupled_2D_ntz(:,i)=nodetrans(coupled_2D_index,2);
% end       sita_part(tt,:,i)

%3D待优化区内节点索引匹配
trshd=0.2;         %比较关键
coupled_3D_index(:,6)=coupled_2D_index;
for i=1:length(coupled_2D_index)
    for j=1:zn
        sita_temp=sita_part2(:,:,j);
        temp2_r_index=find(sita_temp(:,1)<(nodetrans(coupled_2D_index(i),1)+trshd)&(nodetrans(coupled_2D_index(i),1)-trshd)<sita_temp(:,1));
        temp2_z_index=find(sita_temp(:,2)<(nodetrans(coupled_2D_index(i),2)+trshd*2)&(nodetrans(coupled_2D_index(i),2)-trshd*2)<sita_temp(:,2));
        
        for k=1:length(temp2_r_index)   %通用性提高
            if find(temp2_r_index(k)==temp2_z_index(:)) %    && find(temp_rz_index(j)==temp_z_index(:))              
                coupled_3D_index(i,j)=sita_temp(temp2_r_index(k),3);                         
            end
        end  
    end 
end

%3D补充节(过渡区)节点索引匹配
suply=xlsread('nodedata.xlsx','suply','A1:A32'); %从mapping文件的suply_det变量获取，输入excel
coupled_2D_index=[coupled_2D_index;suply];  %2D补充

suply_r=nodetrans(suply,1);
suply_z=nodetrans(suply,2);

ll=length(coupled_3D_index);
for i=1:length(suply)
	tt=1;    
	temp2_r_index=find(nodetrans(:,1)<(suply_r(i)+trshd)&(suply_r(i)-trshd)<nodetrans(:,1));
	temp2_z_index=find(nodetrans(:,2)<(suply_z(i)+trshd)&(suply_z(i)-trshd)<nodetrans(:,2));
        
	for k=1:length(temp2_r_index)   %通用性提高
        if find(temp2_r_index(k)==temp2_z_index(:)) %    && find(temp_rz_index(j)==temp_z_index(:))
            coupled_3D_temp(1,tt)=temp2_r_index(k);               
            tt=tt+1;               
        end
    end
    
    ll=ll+1;     
	coupled_3D_index(ll,:)=coupled_3D_temp;
end

%% 筛选除待更新节点外的固定节点
nodefixed=nodeinitial;  
nodefixed(coupled_3D_index,:)=[]; 

%% 求解chebyshev点对应的Z坐标值
left_R=nodetrans(left_index,1);%left(:,3);
left_Z=nodetrans(left_index,2);%left(:,2);
left_cheby_z=pchip(left_R,left_Z,left_cheby_r);     %R的值不能重复   作为PSO初始边界
% left_cheby_z=interp1(left_R,left_Z,left_cheby_r,'Nearest');

right_R=nodetrans(right_index,1);
right_Z=nodetrans(right_index,2);
right_cheby_z=pchip(right_R,right_Z,right_cheby_r); %R的值不能重复   作为PSO初始边界
% right_cheby_z=interp1(right_R,right_Z,right_cheby_r,'Nearest');

%% lagrange插值
syms t
Fl=lagrange(left_cheby_r,left_cheby_z);         %输出多项式
left_chelag_r=double(vpa(subs(Fl,t,a:b1),6));   %代入求解

Fr=lagrange(right_cheby_r,right_cheby_z);       %输出多项式
right_chelag_r=double(vpa(subs(Fr,t,a:b2),6));  %代入求解

%% PSOA 粒子初始化
left_min_z=left_cheby_z-[0.5;0.6;0.7;0.8;0.8;0.8;0.7;0.6;0.5];      %参数搜索范围
left_max_z=left_cheby_z+[2;3;4.5;5;5;4.5;4;3;2];
right_min_z=right_cheby_z-[2;3;4.5;5;5;5;4.5;3;2];      %参数搜索范围
right_max_z=right_cheby_z+[0.5;0.6;0.7;0.8;0.8;0.8;0.7;0.6;0.5];

MinX=[left_min_z;right_min_z];
MaxX=[left_max_z;right_max_z];

Vmax=1;              %限定速度的范围
Vmin=-1;

Size=50;        %粒子组数
CodeL=length(left_min_z)+length(right_min_z);         %参数个数

c1=1.9;c2=2.0;c3=0.1;c4=0.15;   %学习因子
wmax=0.9;wmin=0.4; %惯性权重范围
G=600;         %最大迭代次数

for i=1:G
    w(i)=wmax-((wmax-wmin)/G)*i;  
end

for i=1:1:CodeL/2 %参数个数
    X(:,i)=left_cheby_z(i)+(left_max_z(i)-left_min_z(i))*rand(Size,1);          
    X(:,CodeL/2+i)=right_cheby_z(i)+(right_max_z(i)-right_min_z(i))*rand(Size,1);                                    
    v(:,i)=Vmin+(Vmax-Vmin)*rand(Size,1);                                     
    v(:,CodeL/2+i)=Vmin+(Vmax-Vmin)*rand(Size,1);
end

for i=1:1:Size
	ssaleft=ssa([nodetrans(2956,2),nodetrans(2947,2),X(i,1:CodeL/2),nodetrans(2623,2),nodetrans(2624,2)]); 
	X(i,1:CodeL/2)=ssaleft(3:length(ssaleft)-2);
	ssaright=ssa([nodetrans(2730,2),nodetrans(2734,2),X(i,CodeL/2+1:CodeL),nodetrans(2681,2),nodetrans(2684,2)]);
	X(i,CodeL/2+1:CodeL)=ssaright(3:length(ssaright)-2);
end

%% 切比雪夫转换为节点坐标
transformation

%% 过渡段网格调整
join=xlsread('nodedata.xlsx','join','A1:D20'); 
join_lu_index=join(1:4,:);
join_ru_index=join(5:8,:);
join_ld_index=join(9:14,:);
join_rd_index=join(15:20,:);

%过渡段矩阵初始化(左上，右上，左下，右下）
join_lu_z=zeros(size(join_lu_index,1),size(join_lu_index,2),Size);
join_ru_z=zeros(size(join_ru_index,1),size(join_ru_index,2),Size);
join_ld_z=zeros(size(join_ld_index,1),size(join_ld_index,2),Size);
join_rd_z=zeros(size(join_rd_index,1),size(join_rd_index,2),Size);

join_lu_r=zeros(size(join_lu_index,1),size(join_lu_index,2),Size);
join_ru_r=zeros(size(join_ru_index,1),size(join_ru_index,2),Size);
join_ld_r=zeros(size(join_ld_index,1),size(join_ld_index,2),Size);
join_rd_r=zeros(size(join_rd_index,1),size(join_rd_index,2),Size);

%初始赋值，根据相对节点，给予原始Z值
for i=1:size(join,2) %按列排列
	join_lu_z(:,i,1)=nodetrans(join_lu_index(:,i),2);
	join_ru_z(:,i,1)=nodetrans(join_ru_index(:,i),2);
end

for i=1:size(join,2) %按列排列
	join_ld_z(:,i,1)=nodetrans(join_ld_index(:,i),2);
	join_rd_z(:,i,1)=nodetrans(join_rd_index(:,i),2);
end

%补齐固定节点（调试过程可注释掉）(也可初始化中计算，不必重复计算)
for i=1:Size
    join_lu_z(1,:,i)=join_lu_z(1,:,1);
    join_ru_z(1,:,i)=join_ru_z(1,:,1);
    join_ld_z(size(join_ld_index,1),:,i)=join_ld_z(size(join_ld_index,1),:,1);
    join_rd_z(size(join_rd_index,1),:,i)=join_rd_z(size(join_rd_index,1),:,1);
    
    join_lu_z(:,1,i)=join_lu_z(:,1,1);
    join_ru_z(:,size(join_ru_index,2),i)=join_ru_z(:,size(join_ru_index,2),1);
    join_ld_z(:,1,i)=join_ld_z(:,1,1);
    join_rd_z(:,size(join_rd_index,2),i)=join_rd_z(:,size(join_rd_index,2),1);   
end

%函数边界的内部节点更新（交点求解,参考线求解）
ref_p=2;

% 给原始r赋值
for i=1:size(join_lu_index,2)
    for j=1:Size %按列排列
        join_lu_r(:,i,j)=nodetrans(join_lu_index(:,i),1);
        join_ru_r(:,i,j)=nodetrans(join_ru_index(:,i),1);
        join_ld_r(:,i,j)=nodetrans(join_ld_index(:,i),1);
        join_rd_r(:,i,j)=nodetrans(join_rd_index(:,i),1);
    end
end

tirt=cell(Size,ref_p);
tizt=cell(Size,ref_p);

regrid

%% 节点映射
clear x1 y1 z1
mapping

%% 节点位置输出
export

%% ANSYS求解
system('"D:\Tools\ANSYS Inc\v202\ansys\bin\winx64\ANSYS202" -b -p ane3fl -i C:\Users\DongZZ\Desktop\optResearch\mixedcode\Main_mix.f -o D:\DZZ.out')   

%% PSOA 最优解集初始化
load matapdlrst.txt         %导入求解结果（FROM APDL）
body_eqv=matapdlrst(:,1);
% root_eqv=matapdlrst(:,2);
disc_volu=matapdlrst(:,2);

for i=1:length(body_eqv)
	qq(i)=max(0,body_eqv(i)-limit_body_eqv);          %罚函数
        
	if qq(i)<10
        score(i)=disc_volu(i)+qq(i)^2;        
	elseif qq(i)>=10 && qq(i)<100
        score(i)=disc_volu(i)+qq(i)^2.5;           
    else
        score(i)=disc_volu(i)+qq(i)^3;   
	end      
end

%  个体最佳初始化
Xhb=X;            
scorehb=score;

[~,tm]=min(score);
BestS=X(tm,:);
scorehgb=scorehb(tm);

hisBestScore=scorehgb;

%% 粒子群更新主程序
BestJ=zeros(G,1);
optPath=zeros(G,CodeL);
% 
disPN=scorehgb-1;    %用于显示结果
% hisBestScore=scorehgb;
hwait=waitbar(0,'总进度>>>>>>>>'); 
for kg=1:1:G             %最大迭代次数
    BestJ(kg)=vpa(scorehgb/hisBestScore,4);
    optPath(kg,:)=BestS;
    for i=1:Size         %粒子组数
        v(i,:)=w(kg)*v(i,:)+c1*rand*(Xhb(i,:)-X(i,:))+c2*rand*(BestS-X(i,:));%+c3*rand*([left_max_z;right_min_z]'-X(i,:));%c3*rand*([(X(i,1:CodeL/2)+X(i,CodeL/2+1:CodeL))./2,(X(i,1:CodeL/2)+X(i,CodeL/2+1:CodeL))./2]-X(i,:));%+c4*rand*([nodetrans(couplednodes(:,1),2);nodetrans(couplednodes(:,mn+2),2)]'-X(i,:));       
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
    end  
    
    if rem(kg,15)==1 && kg<G/2
        for i=1:1:Size
            ssaleft=ssa([nodetrans(2956,2),nodetrans(2947,2),X(i,1:CodeL/2),nodetrans(2623,2),nodetrans(2624,2)]); 
            X(i,1:CodeL/2)=ssaleft(3:length(ssaleft)-2);
            ssaright=ssa([nodetrans(2730,2),nodetrans(2734,2),X(i,CodeL/2+1:CodeL),nodetrans(2681,2),nodetrans(2684,2)]);
            X(i,CodeL/2+1:CodeL)=ssaright(3:length(ssaright)-2);
        end
        save rst.mat
    end   

    %输出粒子位置
    transformation
    regrid
    mapping
    export

    %ANSYS求解
    system('"D:\Tools\ANSYS Inc\v202\ansys\bin\winx64\ANSYS202" -b -p ane3fl -i C:\Users\DongZZ\Desktop\optResearch\mixedcode\Main_mix.f -o D:\DZZ.out')
        
    load matapdlrst.txt         %导入求解结果（FROM APDL）
	body_eqv=matapdlrst(:,1);
	disc_volu=matapdlrst(:,2);
        
    %计算适应度值
    for i=1:length(body_eqv)
        qq(i)=max(0,body_eqv(i)-limit_body_eqv);          %罚函数
        
        if qq(i)<10
            score(i)=disc_volu(i)+qq(i)^2;        
        elseif qq(i)>=10 && qq(i)<100
            score(i)=disc_volu(i)+qq(i)^2.5;           
        else
            score(i)=disc_volu(i)+qq(i)^3;   
        end   
    end
    
	%更新全局、局部最优解
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
end

% delete(gcp('nocreate')) 
BestJ
optPath
