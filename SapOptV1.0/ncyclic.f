!节点循环对称
esel,s,type,,face1   
cm,a_low,node
esel,s,type,,face2 
nsle,all
cm,a_high,node
esel,s,type,,1 	!类型1网格为对称面的体网格
cyclic,36,360/36   !,1,a,0,  !共36个叶片，角度为360/36