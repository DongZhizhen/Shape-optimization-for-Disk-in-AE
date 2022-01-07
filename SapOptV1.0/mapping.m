%% 从边界1D推广到（极坐标）2D
for j=1:left_coup_index(2)-left_coup_index(1)+1
    delt_r(j)=nodetrans(right_index(right_coup_index(1)+j-1),1)-nodetrans(left_index(left_coup_index(1)+j-1),1);    % right_pso_z(right_coup_index(1)+j-1,i)-left_pso_z(left_coup_index(1)+j-1,i);
end

delt_z=zeros(Size,left_coup_index(2)-left_coup_index(1)+1);
for i=1:1:Size 
    for j=1:CodeL/2
        delt_z(i,:)=X(i,CodeL/2+1:CodeL)'-X(i,1:CodeL/2)';
    end
end

for i=1:1:Size 
    for j=1:CodeL/2
        for k=1:mn
            mid_point_r(i,j,k)=nodetrans(left_index(left_coup_index(1)+j-1),1)+k*delt_r(j)/(mn+1);
            mid_point_z(i,j,k)=X(i,j)+k*delt_z(i,j)/(mn+1);
        end
    end
end

update_2D_r=zeros(Size,CodeL/2,mn+2);% R坐标初始化
update_2D_z=zeros(Size,CodeL/2,mn+2);% Z坐标初始化

for i=1:Size
    update_2D_r(i,:,1)=nodetrans(left_index(left_coup_index(1):1:left_coup_index(2)))';   %左边节点
end
update_2D_z(:,:,1)=X(:,1:CodeL/2);   %左边节点

update_2D_r(:,:,2:mn+1)=mid_point_r;   %中间节点
update_2D_z(:,:,2:mn+1)=mid_point_z;   %中间节点

for i=1:Size
    update_2D_r(i,:,mn+2)=nodetrans(right_index(right_coup_index(1):1:right_coup_index(2)))';   %右边节点
end
update_2D_z(:,:,mn+2)=X(:,CodeL/2+1:CodeL);   %右边节点

%% 极坐标2D坐标关联节点（极坐标）2D
sita=unique(nodetrans(:,4));  %获取角度坐标值

single_2D_pol=zeros((left_coup_index(2)-left_coup_index(1)+1)*(mn+2),4);
update_2D_pol=zeros(Size,(left_coup_index(2)-left_coup_index(1)+1)*(mn+2),4);

for i=1:Size
    temp=1;
	for k=1:mn+2
        for j=1:left_coup_index(2)-left_coup_index(1)+1   %单粒子节点坐标跟新
            single_2D_pol(temp,:)=[couplednodes(j,k),sita(1),update_2D_r(i,j,k),update_2D_z(i,j,k)];
            temp=temp+1;
        end
	end
    
    for j=1:(left_coup_index(2)-left_coup_index(1)+1)*(mn+2)  %全部粒子节点坐标跟新
        for k=1:4
            update_2D_pol(i,j,k)=single_2D_pol(j,k);
        end
    end
end

%% 从2D推广到（极坐标）3D
sita=unique(nodetrans(:,4));

dl_sita=(max(nodetrans(:,4))-min(nodetrans(:,4)))/zn;

single_3D_pol=zeros((left_coup_index(2)-left_coup_index(1)+1)*(mn+2)*(zn+1),4);
update_3D_pol=zeros(Size,(left_coup_index(2)-left_coup_index(1)+1)*(mn+2)*(zn+1),4);

for i=1:Size
    temp=1;    
    tempd=1;    
    for j=1:left_coup_index(2)-left_coup_index(1)+1   %单粒子节点坐标跟新
        for k=1:mn+2
            for d=1:zn+1
                single_3D_pol(tempd,:)=[coupled_3D(temp,d),sita(1)+(d-1)*dl_sita,update_2D_pol(i,temp,3),update_2D_pol(i,temp,4)];    
                tempd=tempd+1;
            end
            temp=temp+1;
        end
    end    
    
    for j=1:(left_coup_index(2)-left_coup_index(1)+1)*(mn+2)*(zn+1)  %全部粒子节点坐标跟新
        for k=1:4
            update_3D_pol(i,j,k)=single_3D_pol(j,k);
        end
    end
end

%% （极坐标）3D转换为（笛卡尔）3D
update_3D_cart=zeros(Size,(left_coup_index(2)-left_coup_index(1)+1)*(mn+2)*(zn+1),4);

clear x y z
for i=1:Size
    [x,y,z] = pol2cart(update_3D_pol(i,:,2),update_3D_pol(i,:,3),update_3D_pol(i,:,4)); 
    update_temp=[update_3D_pol(i,:,1)',x',z',y'];    %对应模型坐标（上节导过来）
    
    for j=1:(left_coup_index(2)-left_coup_index(1)+1)*(mn+2)*(zn+1)
        for k=1:4
            update_3D_cart(i,j,k)=update_temp(j,k);
        end
    end
end