program domain_decomposition
real xmin,ymin,zmin,xmax,ymax,zmax,aux
real xmin2,ymin2,zmin2,xmax2,ymax2,zmax2,rmax
real, allocatable :: my_array(:,:),wgt1(:)
integer, allocatable :: gal(:),region(:)
integer Ndata,i,j,nregion,buffer
character*100 c 


CALL getarg(1, c)
read(c,*) rmax
print*, rmax

Ndata=1
open(11,file='allpoints',status='unknown')
21 continue
read(11,*,err=21,end=22) aux
Ndata=Ndata+1
goto 21
22 close(11)

allocate(my_array(3,Ndata))
allocate(wgt1(Ndata))
allocate(gal(Ndata))
allocate(region(Ndata))

i=1
open(11,file='allpoints',status='unknown')
25 continue
read(11,*,err=25,end=26)my_array(1:3,i),wgt1(i),gal(i),aux,region(i)
i=i+1
goto 25
26 close(11)

nregion=maxval(region)

print*, "creating ",nregion," regions"

do j=1,nregion

xmin=minval(my_array(1,:),mask=region.eq.j)
ymin=minval(my_array(2,:),mask=region.eq.j)
zmin=minval(my_array(3,:),mask=region.eq.j)
xmax=maxval(my_array(1,:),mask=region.eq.j)
ymax=maxval(my_array(2,:),mask=region.eq.j)
zmax=maxval(my_array(3,:),mask=region.eq.j)

print*,xmin,xmax,ymin,ymax,zmin,zmax
xmin2=xmin-rmax
ymin2=ymin-rmax
zmin2=zmin-rmax
xmax2=xmax+rmax
ymax2=ymax+rmax
zmax2=zmax+rmax
print*,xmin2,xmax2,ymin2,ymax2,zmin2,zmax2

write (c, "(I3,A10)") j,".loadnodes"
c=adjustl(c)
print*, c!//"test.node"

open(unit=11, file=trim(c),status='unknown')
do i=1,Ndata

if (my_array(1,i)>xmin2 .and. my_array(1,i)<xmax2 .and. &
&   my_array(2,i)>ymin2 .and. my_array(2,i)<ymax2 .and. &
&   my_array(3,i)>zmin2 .and. my_array(3,i)<zmax2) then


buffer=1
if (my_array(1,i)>=xmin .and. my_array(1,i)<=xmax .and. &
& my_array(2,i)>=ymin .and. my_array(2,i)<=ymax .and. &
& my_array(3,i)>=zmin .and. my_array(3,i)<=zmax) buffer=0

write(11,*) my_array(:,i),wgt1(i),gal(i),buffer

endif
enddo
close(11)

enddo

deallocate(my_array,wgt1,gal,region)

end program domain_decomposition
