!��������
/clear
size=50         !���Ӹ���
objs=2        !�Ż�Ŀ���������Ŀ�꣩

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
	nread,chebyschemes_%i%,txt      !�������ڵ㣬matlab������
	et,1,185       !���嵥Ԫ����
	et,2,154

	*if,i,lt,1.5,then,
		eread,elemcal,elem      !�������ڵ㣬matlab������
		*use,materials.f
		*use,ncyclic.f
	*endif

	nplot
	eplot 

	/solu
	*use,solve.f

	!���ӻ�
	/post1
	set,first
	/dscale,all,1.0 
	/EFACET,1   
	plesol,s,eqv, 0,1.0

	*use,nprobe.f

*enddo

*use,writing.f     !������







