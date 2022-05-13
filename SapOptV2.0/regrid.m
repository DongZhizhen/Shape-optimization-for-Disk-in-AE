%角节点更新
join_lu_z(end,end,:)=left_border_z(1,:);
join_ru_z(end,1,:)=right_border_z(1,:);
join_ld_z(1,end,:)=left_border_z(end,:);
join_rd_z(1,1,:)=right_border_z(end,:);

%四直线边的内部节点更新
zcu=3;%纵向3层网格
for i=1:zcu-1
    for j=1:Size  
        join_lu_z(i+1,end,j)=join_lu_z(1,end,1)-i*(join_lu_z(1,end,1)-join_lu_z(end,end,j))/zcu;
        join_ru_z(i+1,1,j)=join_ru_z(1,1,1)-i*(join_ru_z(1,1,1)-join_ru_z(end,1,j))/zcu;    
    end
end

zcd=5;%纵向5层网格
for i=1:zcd-1
    for j=1:Size  
%         join_ld_z(i+1,size(join_ld_index,2),j)=join_ld_z(1,size(join_ld_index,2),1)-i*(join_ld_z(1,size(join_ld_index,2),1)-join_ld_z(size(join_ld_index,1),size(join_ld_index,2),j))/zcd;
%         join_rd_z(i+1,1,j)=join_rd_z(1,1,j)-i*(join_rd_z(1,1,j)-join_rd_z(size(join_rd_index,1),1,1))/zcd;
        join_ld_z(i+1,end,j)=join_ld_z(1,end,j)-i*(join_ld_z(1,end,j)-join_ld_z(end,end,1))/zcd;
        join_rd_z(i+1,1,j)=join_rd_z(1,1,j)-i*(join_rd_z(1,1,j)-join_rd_z(end,1,1))/zcd;
    end
end

%函数边界的内部节点更新（交点求解,参考点求解）
for i=1:size(join,2)-2
    for j=1:Size  
        join_lu_z(size(join_lu_index,1),i+1,j)=join_lu_z(size(join_lu_index,1),1,1)+i*(join_lu_z(size(join_lu_index,1),size(join_lu_index,2),j)-join_lu_z(size(join_lu_index,1),1,1))/(size(join,2)-1);
        join_ru_z(size(join_ru_index,1),i+1,j)=join_ru_z(size(join_ru_index,1),1,j)+i*(join_ru_z(size(join_ru_index,1),size(join_ru_index,2),1)-join_ru_z(size(join_ru_index,1),1,j))/(size(join,2)-1);
    end
end

for i=1:size(join,2)-2
    for j=1:Size  
        join_ld_z(1,i+1,j)=join_ld_z(1,1,1)+i*(join_ld_z(1,size(join_ld_index,2),j)-join_ld_z(1,1,1))/(size(join,2)-1);                           
        join_rd_z(1,i+1,j)=join_rd_z(1,1,j)+i*(join_rd_z(1,size(join_rd_index,2),1)-join_rd_z(1,1,j))/(size(join,2)-1); 
    end
end

% %函数边界的内部节点更新（交点求解,参考线求解）
%左上1,2
ref_x(:,1)=[nodetrans(join_lu_index(1,ref_p),1);nodetrans(join_lu_index(end,ref_p),1)];
ref_x(:,ref_p)=[nodetrans(join_lu_index(1,ref_p+1),1);nodetrans(join_lu_index(end,ref_p+1),1)];

for i=1:Size
    ref_y(:,1,i)=[join_lu_z(1,ref_p,i);join_lu_z(end,ref_p,i)];    
    ref_y(:,ref_p,i)=[join_lu_z(1,ref_p+1,i);join_lu_z(end,ref_p+1,i)];
    
    for j=1:ref_p             %横坐标1~2为方程1和方程2
        ref_lu_k(j,i)=(ref_y(1,j,i)-ref_y(2,j,i))/(ref_x(1,j)-ref_x(2,j));
        ref_lu_b(j,i)=ref_y(1,j,i)-ref_lu_k(j,i)*ref_x(1,j);
    end
end

%右上1,2
ref_x(:,1)=[nodetrans(join_ru_index(1,ref_p),1);nodetrans(join_ru_index(end,ref_p),1)];
ref_x(:,ref_p)=[nodetrans(join_ru_index(1,ref_p+1),1);nodetrans(join_ru_index(end,ref_p+1),1)];

for i=1:Size
    ref_y(:,1,i)=[join_ru_z(1,ref_p,i);join_ru_z(end,ref_p,i)];    
    ref_y(:,ref_p,i)=[join_ru_z(1,ref_p+1,i);join_ru_z(end,ref_p+1,i)];
    for j=1:ref_p             %横坐标1~2为方程1和方程2
        ref_ru_k(j,i)=(ref_y(1,j,i)-ref_y(2,j,i))/(ref_x(1,j)-ref_x(2,j));
        ref_ru_b(j,i)=ref_y(1,j,i)-ref_ru_k(j,i)*ref_x(1,j);
    end
end

%左下1,2
ref_x(:,1)=[nodetrans(join_ld_index(1,ref_p),1);nodetrans(join_ld_index(end,ref_p),1)];
ref_x(:,ref_p)=[nodetrans(join_ld_index(1,ref_p+1),1);nodetrans(join_ld_index(end,ref_p+1),1)];

for i=1:Size
    ref_y(:,1,i)=[join_ld_z(1,ref_p,i);join_ld_z(end,ref_p,i)];    
    ref_y(:,ref_p,i)=[join_ld_z(1,ref_p+1,i);join_ld_z(end,ref_p+1,i)];
    for j=1:ref_p             %横坐标1~2为方程1和方程2
        ref_ld_k(j,i)=(ref_y(1,j,i)-ref_y(2,j,i))/(ref_x(1,j)-ref_x(2,j));
        ref_ld_b(j,i)=ref_y(1,j,i)-ref_ld_k(j,i)*ref_x(1,j);
    end
end

%右下1,2
ref_x(:,1)=[nodetrans(join_rd_index(1,ref_p),1);nodetrans(join_rd_index(end,ref_p),1)];
ref_x(:,ref_p)=[nodetrans(join_rd_index(1,ref_p+1),1);nodetrans(join_rd_index(end,ref_p+1),1)];

for i=1:Size
    ref_y(:,1,i)=[join_rd_z(1,ref_p,i);join_rd_z(end,ref_p,i)];    
    ref_y(:,ref_p,i)=[join_rd_z(1,ref_p+1,i);join_rd_z(end,ref_p+1,i)];
    for j=1:ref_p             %横坐标1~2为方程1和方程2
        ref_rd_k(j,i)=(ref_y(1,j,i)-ref_y(2,j,i))/(ref_x(1,j)-ref_x(2,j));
        ref_rd_b(j,i)=ref_y(1,j,i)-ref_rd_k(j,i)*ref_x(1,j);
    end
end

% %函数边界的内部节点更新（交点求解,节点更新）
syms iz ir
%左上1,2
parfor i=1:Size
    for j=1:ref_p
        S=vpasolve(iz==subs(Fl_pso(i),t,ir),iz==ref_lu_k(j,i)*ir+ref_lu_b(j,i),iz,ir);
        tirt{i,j}=S.ir(imag(S.ir)==0&real(S.iz)<=0);
        tizt{i,j}=S.iz(imag(S.iz)==0&real(S.iz)<=0);
    end
end

for i=1:Size
    for j=1:ref_p        
        try         %尝试求交点
            if length(tirt{i,j})==1
                join_lu_r(end,j+1,i)=tirt{i,j};
                join_lu_z(end,j+1,i)=tizt{i,j};
            else
                tir=tirt{i,j};
                tiz=tizt{i,j};
                join_lu_r(end,j+1,i)=tir(tir>(b1-6)&tir<b1);
                join_lu_z(end,j+1,i)=tiz(tir>(b1-6)&tir<b1);
            end
        catch       %报错采用折衷方案，参考点系数0.7，切比雪夫函数0.3
            join_lu_r(end,j+1,i)=nodetrans(join_lu_index(end,j+1),1);
            join_lu_z(end,j+1,i)=0.3*double(vpa(subs(Fl_pso(i),t,nodetrans(join_lu_index(end,j+1),1)),6))+0.7*join_lu_z(end,j+1,i);
        end
        
        for k=2:size(join_lu_index,1)-1
            join_lu_r(k,j+1,i)=join_lu_r(1,j+1,i)-(k-1)*(join_lu_r(1,j+1,i)-join_lu_r(end,j+1,i))/(size(join_lu_index,1)-1);
        end
    end
end

%右上1,2
parfor i=1:Size
    for j=1:ref_p
        S=vpasolve(iz==subs(Fr_pso(i),t,ir),iz==ref_ru_k(j,i)*ir+ref_ru_b(j,i),iz,ir);
        tirt{i,j}=S.ir(imag(S.ir)==0&real(S.iz)>=0);
        tizt{i,j}=S.iz(imag(S.iz)==0&real(S.iz)>=0);
    end
end

for i=1:Size
    for j=1:ref_p        
        try         %尝试求交点
            if length(tirt{i,j})==1
                join_ru_r(end,j+1,i)=tirt{i,j};
                join_ru_z(end,j+1,i)=tizt{i,j};
            else
                tir=tirt{i,j};
                tiz=tizt{i,j}; 
                join_ru_r(end,j+1,i)=tir(tir>(b2-6)&tir<b2);
                join_ru_z(end,j+1,i)=tiz(tir>(b2-6)&tir<b2);   
            end
        catch       %报错采用折衷方案，参考点系数0.7，切比雪夫函数0.3    
            join_ru_r(end,j+1,i)=nodetrans(join_ru_index(end,j+1),1);
            join_ru_z(end,j+1,i)=0.3*double(vpa(subs(Fr_pso(i),t,nodetrans(join_ru_index(end,j+1),1)),6))+0.7*join_ru_z(end,j+1,i);        
        end
        
        for k=2:size(join_ru_index,1)-1
            join_ru_r(k,j+1,i)=join_ru_r(1,j+1,i)-(k-1)*(join_ru_r(1,j+1,i)-join_ru_r(end,j+1,i))/(size(join_ru_index,1)-1);
        end
    end
end

%左下1,2
parfor i=1:Size
    for j=1:ref_p
        S=vpasolve(iz==subs(Fl_pso(i),t,ir),iz==ref_ld_k(j,i)*ir+ref_ld_b(j,i),iz,ir);
        tirt{i,j}=S.ir(imag(S.ir)==0&real(S.iz)<=0);
        tizt{i,j}=S.iz(imag(S.iz)==0&real(S.iz)<=0);
    end
end

for i=1:Size
    for j=1:ref_p
        try         %尝试求交点
            if length(tirt{i,j})==1
                join_ld_r(1,j+1,i)=tirt{i,j};
                join_ld_z(1,j+1,i)=tizt{i,j};
            else
                tir=tirt{i,j};
                tiz=tizt{i,j}; 
                join_ld_r(1,j+1,i)=tir(tir<(a+6)&tir>a);
                join_ld_z(1,j+1,i)=tiz(tir<(a+6)&tir>a);           
            end
        catch       %报错采用折衷方案，参考点系数0.7，切比雪夫函数0.3   
            join_ld_r(1,j+1,i)=nodetrans(join_ld_index(1,j+1),1);
            join_ld_z(1,j+1,i)=0.3*double(vpa(subs(Fl_pso(i),t,nodetrans(join_ld_index(1,j+1),1)),6))+0.7*join_ld_z(1,j+1,i); 
        end
                            
        for k=2:size(join_ld_index,1)-1
            join_ld_r(k,j+1,i)=join_ld_r(1,j+1,i)-(k-1)*(join_ld_r(1,j+1,i)-join_ld_r(end,j+1,i))/(size(join_ld_index,1)-1);
        end
    end
end

%右下1,2
parfor i=1:Size
    for j=1:ref_p
        S=vpasolve(iz==subs(Fr_pso(i),t,ir),iz==ref_rd_k(j,i)*ir+ref_rd_b(j,i),iz,ir);
        tirt{i,j}=S.ir(imag(S.ir)==0&real(S.iz)<=0);
        tizt{i,j}=S.iz(imag(S.iz)==0&real(S.iz)<=0);
    end
end

for i=1:Size
    for j=1:ref_p        
        try         %尝试求交点
            if length(tirt{i,j})==1
                join_rd_r(1,j+1,i)=tirt{i,j};
                join_rd_z(1,j+1,i)=tizt{i,j};
            else
                tir=tirt{i,j};
                tiz=tizt{i,j}; 
                join_rd_r(1,j+1,i)=tir(tir<(a+6)&tir>a);
                join_rd_z(1,j+1,i)=tiz(tir<(a+6)&tir>a);    
            end
        catch       %报错采用折衷方案，参考点系数0.7，切比雪夫函数0.3  
            join_rd_r(1,j+1,i)=nodetrans(join_rd_index(1,j+1),1);
            join_rd_z(1,j+1,i)=0.3*double(vpa(subs(Fr_pso(i),t,nodetrans(join_rd_index(1,j+1),1)),6))+0.7*join_rd_z(1,j+1,i);    
        end   
                    
        for k=2:size(join_rd_index,1)-1
            join_rd_r(k,j+1,i)=join_rd_r(1,j+1,i)-(k-1)*(join_rd_r(1,j+1,i)-join_rd_r(end,j+1,i))/(size(join_rd_index,1)-1);
        end
    end
end
delete(gcp('nocreate')) 

%过渡段中间节点坐标计算
%上层网格zcu共享
for k=1:size(join,2)-2
    for i=1:zcu-1
        for j=1:Size  
            join_lu_z(i+1,k+1,j)=join_lu_z(1,k+1,1)-i*(join_lu_z(1,k+1,1)-join_lu_z(size(join_lu_index,1),k+1,j))/zcu;
            join_ru_z(i+1,k+1,j)=join_ru_z(1,k+1,1)-i*(join_ru_z(1,k+1,1)-join_ru_z(size(join_ru_index,1),k+1,j))/zcu;    
        end
    end
end
%下层网格zcd共享
for k=1:size(join,2)-2
    for i=1:zcd-1
        for j=1:Size  
            join_ld_z(i+1,k+1,j)=join_ld_z(1,k+1,j)-i*(join_ld_z(1,k+1,j)-join_ld_z(size(join_ld_index,1),k+1,1))/zcd;
            join_rd_z(i+1,k+1,j)=join_rd_z(1,k+1,j)-i*(join_rd_z(1,k+1,j)-join_rd_z(size(join_ld_index,1),k+1,1))/zcd;
        end
    end
end
