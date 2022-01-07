!迭代计算
/clear
size=50         !粒子个数
objs=2        !优化目标数（输出目标）

/prep7
/nerr,2,50000

csys,5

/prep7
nsel,all
nrotat,all

*dim,disc_para,array,size,objs

*do,i,1,size
!/clear
!i=2
	/prep7
	nread,chebyschemes_%i%,txt      !读入计算节点，matlab计算结果
	et,1,185       !定义单元类型
	et,2,154

	*if,i,lt,1.5,then,
		eread,elemcal,elem      !读入计算节点，matlab计算结果
		*use,materials.f
		*use,ncyclic.f
	*endif

	nplot
	eplot 

	/solu
	*use,solve.f

	!可视化
	/post1
	set,first
	/dscale,all,1.0 
	/EFACET,1   
	plesol,s,eqv, 0,1.0

	*use,nprobe.f

*enddo

*use,writing.f     !输出结果







