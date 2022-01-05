!*****************************************************************************
!
!  Main_07.f
!
!  Description:
!      This program is used to run the ANSYS software. 
!      The purpose is to load the grid, calculation the model and output the results.
!
!  All the codes had been written by Zhizhen Dong in 2020-2021
!
!****************************************************************************

!��������
/clear
size=50         !���Ӹ���
objs=3        !�Ż�Ŀ���������Ŀ�꣩

/prep7
/nerr,2,50000

csys,5
nsel,all
nrotat,all

*dim,disc_para,array,size,objs

*do,i,1,size

/prep7
nread,chebyschemes_%i%,f      !�������ڵ㣬matlab������
et,1,185       !���嵥Ԫ����
et,2,154

*if,i,lt,1.5,then,
	eread,elemcal,elem      !�������ڵ㣬matlab������
	*use,materials.f
	*use,ncyclic.f
*endif

nplot
eplot 

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







