program gramsci
  use kdtree2_module
  use kdtree2_precision_module
  implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    INCLUDE "omp_lib.h"
    
  use node_module, only: node
  use sorting_module
  type(node), dimension(:), allocatable :: output
  
  real(kdkind), dimension(:,:), allocatable :: my_array
  real(kdkind), allocatable ::  v1(:),v2(:),v3(:),v4(:)

  type(kdtree2), pointer :: tree
  real(kdkind), allocatable :: Zddd(:,:,:,:,:),Zddr(:,:,:,:,:),Zdrr(:,:,:,:,:),Zrrr(:,:,:,:,:),crap(:,:,:,:,:)
  real(kdkind), allocatable :: Zdd(:,:),Zdr(:,:),Zrr(:,:),crap2(:,:), wgt1(:)
  integer, allocatable :: gal(:), buffer(:)
  real :: avgal,avran
  real(kdkind) :: ndd,nrr
  integer :: k, i, j, l,m,n,d,chunk,nmu,nbins,iloop
  integer*1 :: mu1,mu2,ind1,ind2,ind3,ind4
  character :: file1*2000,file2*2000, ranfile*2000, outfile*2000
  type(kdtree2_result),allocatable :: resultsb(:)
  integer   ::  Ndata,Nrand,nn1,nn2,nnode,nodeid, thread, threads
  real(kdkind)  :: aux,odr,odtheta,theta,start,finish
  real(kdkind) :: rmin, rmax, dr, vec_dist
  logical :: wgt, RSD, loadran,saveran,cut,DOMPI
  integer :: myid , ntasks , ierr
  real(kdkind), parameter :: pival=3.14159265358979323846d0
#ifdef MPI 
  integer , dimension( MPI_STATUS_SIZE ) :: status
#endif
  integer , parameter :: master=0,msgtag1=11,msgtag2=12,msgtag3=13,msgtag4=14


! ! I n i t i a l i z a t i o n of the MPI environment
myid=0
ntasks=1
thread=0
threads=1

#ifdef MPI 
 call MPI_INIT( ierr )
 call MPI_COMM_RANK( MPI_COMM_WORLD , myid, ierr )
 call MPI_COMM_SIZE( MPI_COMM_WORLD , ntasks , ierr )
#endif

 !print*, 'node ',myid, 'of ',ntasks,' checking in...'

!--- set some default parameters ---
  call default_params()

  if(myid==master) print*,'reading options'
  call parseOptions()

  if(myid==master)  print*,'allocated bins'

  if(cut) then
    call count_files()
  else
    call count_files_2()
  endif
 
  Allocate(my_array(d,Ndata))
  Allocate(wgt1(Ndata))

  allocate(gal(Ndata))
  allocate(buffer(Ndata))
  allocate(v1(d))
  allocate(v2(d))
  allocate(v3(d))
  
  if(cut) then
    call read_files()
  else
    call read_files_2()
  endif


if(myid==0) print*,'building tree '
tree => kdtree2_create(my_array,sort=.true.,rearrange=.true.)     ! this is how you create a tree.
!print*,' built tree on node', myid

if(myid==0) print*,'allocating arrays '

call allocate_arrays()



!!!!!!!!!!!!!! sampling data !!!!!!!!!!!!!!!
allocate(v4(100))
do i=1,100
  call random_number(dr)
  j = floor( dr*Ndata )+1
  !print*,dr,j
  nn2=kdtree2_r_count_around_point(tp=tree,idxin=j,correltime=-1,r2=rmax*rmax)
  v4(i)=nn2*1.d0
  !print*,j,nn2
  dr = SUM(v4)/SIZE(v4)
enddo
!print*,maxval(v4),dr,sqrt(SUM((v4-dr)**2)/SIZE(v4))
!dr=dr+2*sqrt(SUM((v4-dr)**2)/SIZE(v4))
!nn2=floor(dr)+1

!print*, 'size of results:', nn2
print*, 'est. memory usage (Gb): ', Ndata*0.004d0*SUM(v4)/(SIZE(v4)*468911.d0)
print*, 'if this is larger than the RAM in your computer, you should probably kill me'

!!!!!!!!!!!! sampling data !!!!!!!!!!!!!!!
deallocate(v4)
allocate(v4(d))

if(myid==0) print*,'building relationships between nodes'

dr=(rmax)/nbins
odr=1./dr
odtheta=0.5*nmu


allocate(output(Ndata))

!$OMP PARALLEL
!$ threads=OMP_GET_NUM_THREADS()
!$ thread=OMP_GET_THREAD_NUM()
!$OMP END PARALLEL

if(myid==master) print*,'Code running with ',ntasks,' MPI processes each with ',threads,' OMP threads'

call cpu_time(start)
call create_graph()
call cpu_time(finish)
print '("Creating the graph took ",f12.3," seconds.")',(finish-start)/threads


call kdtree2_destroy(tree) 

deallocate(my_array)
deallocate(resultsb) 

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

if(myid==0) print*,'finished building node relationships '

call normalise_counts()
call cpu_time(start)
call query_graph_sorted()
!ranfile=trim(outfile)//".4pcf"
!open(11,file=ranfile,status='unknown')
!do iloop=1,nbins
!  ind1=iloop
!  ind2=iloop
!  ind3=iloop
!  ind4=iloop
!  call query_graph_single4(ind1,ind2,ind3,ind4)
!  Zddd=Zddd/(ndd*(ndd-avgal)*(ndd-2*avgal)*(ndd-3*avgal))
!  Zrrr=Zrrr/(nrr*(nrr-avran)*(nrr-2*avran)*(nrr-3*avran))
!  write(11,'(2(e14.7,1x))') (iloop-1)*dr + 0.5*dr,(Zddd(iloop,iloop,iloop,1,1)-Zrrr(iloop,iloop,iloop,1,1))/Zrrr(iloop,iloop,iloop,1,1)
!enddo
!call cpu_time(finish)
!close(11)
!print '("Querying graph took ",f12.3," seconds.")',(finish-start)/threads

!Zddd=0.0
!Zrrr=0.0

!call cpu_time(start)
!ranfile=trim(outfile)//".3pcf"
!open(11,file=ranfile,status='unknown')
!do iloop=1,10
!  ind1=iloop
!  ind2=iloop
!  ind3=iloop
!  call query_graph_single3(ind1,ind2,ind3)
!  Zddd=Zddd/(ndd*(ndd-avgal)*(ndd-2*avgal))
!  Zrrr=Zrrr/(nrr*(nrr-avran)*(nrr-2*avran))
!  write(11,'(2(e14.7,1x))') (iloop-1)*dr + 0.5*dr, (Zddd(iloop,iloop,iloop,1,1)-Zrrr(iloop,iloop,iloop,1,1))/Zrrr(iloop,iloop,iloop,1,1)
!enddo

call cpu_time(finish)
close(11)
print '("Querying graph took ",f12.3," seconds.")',(finish-start)/threads

if(myid==0) print*,'finished querying the graph'

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

if(myid==0) print*,'collecting results from MPI tasks'

#ifdef MPI
call mpi_collect()
#endif

if(myid==0) print*,'results collected' 

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

if(myid==master) then

if(myid==0) print*,'normalising counts' 

if(cut) then
  continue
else
  call normalise_counts()
endif

if(myid==0) print*,'counts normalised' 

ranfile=trim(outfile)//".3pcf"
open(11,file=ranfile,status='unknown')
do l=1,nbins
    do m=1,nbins
        do n=1,nbins
           do j=1,nmu
             do k=1,nmu

             write(11,'(15(e14.7,1x))')(l-1)*dr,l*dr,(m-1)*dr,m*dr,(n-1)*dr,n*dr,&
             ((float(j)-1.)/odtheta)-1.,(float(j)/odtheta)-1.,((float(k)-1.)/odtheta)-1.,(float(k)/odtheta)-1.,&
             Zddd(l,m,n,j,k),Zddr(l,m,n,j,k),Zdrr(l,m,n,j,k),Zrrr(l,m,n,j,k),&
             (Zddd(l,m,n,j,k)-3*Zddr(l,m,n,j,k)+3*Zdrr(l,m,n,j,k)-Zrrr(l,m,n,j,k))/Zrrr(l,m,n,j,k)

            enddo
           enddo
        enddo
    enddo
enddo
close(11)


ranfile=trim(outfile)//".2pcf"
open(11,file=ranfile,status='unknown')
do l=1,nbins
    do k=1,nmu   
             
        write(11,'(8(e14.7,1x))')(l-1)*dr,l*dr,((float(k)-1.)/odtheta)-1.,(float(k)/odtheta)-1.,&
             Zdd(l,k),Zdr(l,k),Zrr(l,k), (Zdd(l,k)-2*Zdr(l,k)+Zrr(l,k))/Zrr(l,k)

    enddo
enddo
close(11)

endif

! deallocate the memory for the data.


#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

call deallocate_arrays()

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE( ierr )
#endif

deallocate(gal)
deallocate(buffer)

if(myid==master) then
  print*, "Exit... stage left!"
  stop
else
  stop
endif

contains

subroutine create_graph()
real :: start, finish
call cpu_time(start)
!$OMP PARALLEL DO schedule(dynamic)  private(j,vec_dist,v1,v3,v4,nn1,nn2,nnode,v2,i,resultsb,theta) & 
!$OMP& shared(tree,Ndata,myid, output)
do i=1,999

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=rmax*rmax,&
&       nfound=nn2,nalloc=(Ndata),results=resultsb)


!!!!!!NEW!!!!!!
    nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=rmin*rmin)  
    nnode=nn2-nn1!-1
!!!!!!NEW!!!!!!

!!!!!!NEW!!!!!!
!    output(i)%nn=nn2-1
     output(i)%nn=nnode
!!!!!!NEW!!!!!!

!!!!!!NEW!!!!!!
    allocate(output(i)%id(nnode))
    allocate(output(i)%dist(nnode))
    allocate(output(i)%mu(nnode))

!    allocate(output(i)%id(nn2-1))
!    allocate(output(i)%dist(nn2-1))
!    allocate(output(i)%mu(nn2-1))
!!!!!!NEW!!!!!!

    v2=my_array(:,i)

!!!!!!NEW!!!!!!
!   do j=2,nn2,1
    do j=nn1+1,nn2,1
!!!!!!NEW!!!!!!
        v1=my_array(:,resultsb(j)%idx)

        call dist(v1,my_array(:,i),vec_dist)

        if(RSD) then
          v3=0.5*[v2(1)+v1(1),v2(2)+v1(2),v2(3)+v1(3)]
          v4=[v2(1)-v1(1),v2(2)-v1(2),v2(3)-v1(3)] ! connecting vector
          theta=v3(1)*v4(1) + v3(2)*v4(2) + v3(3)*v4(3)
          theta=(theta/((sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)) * (sqrt(v4(1)**2 +v4(2)**2 +v4(3)**2))))
          output(i)%mu(j-nn1)=floor((theta+1.)*odtheta)+1
        else
          output(i)%mu(j-nn1)=1
        endif
!!!!!!NEW!!!!!!
        output(i)%dist(j-nn1)=floor(vec_dist*odr)+1
        !output(i)%mu(j-1)=floor(theta*odtheta)+1
        !output(i)%mu(j-nn1)=floor((theta+1.)*odtheta)+1
        output(i)%id(j-nn1)=resultsb(j)%idx

!        output(i)%dist(j-1)=floor(vec_dist*odr)+1
!        !output(i)%mu(j-1)=floor(theta*odtheta)+1
!        output(i)%mu(j-1)=floor((theta+1.)*odtheta)+1
!        output(i)%id(j-1)=resultsb(j)%idx
!!!!!!NEW!!!!!!
    enddo

!!!!!!NEW!!!!!!
!    call sort2(output(i)%id,output(i)%dist,output(i)%mu,nnode)
    call sort2(output(i)%id,output(i)%dist,output(i)%mu,nn2-1)
!!!!!!NEW!!!!!!
    nn2=0

enddo
!$OMP END PARALLEL DO

call cpu_time(finish)
print '("Creating graph will take ~ ",f10.3," minutes.")',(finish-start)*Ndata/(60.*1000.*threads)
print*, 'If this takes longer than time to drink a coffee, maybe you should give me more CPUs'
print*, 'Or consider decomposing the domain decomposition option with larger N.'

!$OMP PARALLEL DO schedule(dynamic)  private(j,vec_dist,v1,v3,v4,nn1,nn2,nnode,v2,i,resultsb,theta) & 
!$OMP& shared(tree,Ndata,myid, output)
do i=1000,Ndata,1

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=rmax*rmax,&
&       nfound=nn2,nalloc=(Ndata),results=resultsb)

!!!!!!NEW!!!!!!
    nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=rmin*rmin)  
    nnode=nn2-nn1
!!!!!!NEW!!!!!!

!!!!!!NEW!!!!!!
    output(i)%nn=nnode

    allocate(output(i)%id(nnode))
    allocate(output(i)%dist(nnode))
    allocate(output(i)%mu(nnode))

!!!!!!NEW!!!!!!
    
    v2=my_array(:,i)

!!!!!!NEW!!!!!!
    do j=nn1+1,nn2,1
!!!!!!NEW!!!!!!
        v1=my_array(:,resultsb(j)%idx)

        call dist(v1,my_array(:,i),vec_dist)

        if(RSD) then
          v3=0.5*[v2(1)+v1(1),v2(2)+v1(2),v2(3)+v1(3)]
          v4=[v2(1)-v1(1),v2(2)-v1(2),v2(3)-v1(3)] ! connecting vector
          theta=v3(1)*v4(1) + v3(2)*v4(2) + v3(3)*v4(3)
          theta=(theta/((sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)) * (sqrt(v4(1)**2 +v4(2)**2 +v4(3)**2))))
          output(i)%mu(j-nn1)=floor((theta+1.)*odtheta)+1
        else
          output(i)%mu(j-nn1)=1
        endif

        output(i)%dist(j-nn1)=floor(vec_dist*odr)+1
        output(i)%id(j-nn1)=resultsb(j)%idx

    enddo
    
!!!!!!NEW!!!!!!
    call sort2(output(i)%id,output(i)%dist,output(i)%mu,nnode)
!!!!!!NEW!!!!!!

    nn2=0

enddo
!$OMP END PARALLEL DO

end subroutine


subroutine query_graph_single4 (ind1,ind2,ind3,ind4)
integer*1 :: ind1,ind2,ind3,ind4
integer :: count, m1,m2,m3,m4,id2,id3,id4,m

mu1=1
mu2=1


!$OMP PARALLEL DO schedule(dynamic)  private(i,j,k,l,m,m1,m2,m3,m4,id2,id3,id4,count) & ! , 
!$OMP& shared(wgt1,gal,output,buffer,mu1,mu2,ind1,ind2,ind3,ind4)&
!$OMP& reduction(+:Zddd,Zddr,Zdrr,Zrrr)
do i=1,Ndata,1             ! begin loop over all data (and random) point
  if(buffer(i)==1) cycle ! but skip if in buffer region

    m1=output(i)%nn       ! open node A, m1 = # of neighbors

    do j=1,m1             ! loop over neighbor list

      if(output(i)%dist(j).eq.ind1) then !if 1st bin satisfied
        id2=output(i)%id(j)
        m2=output(id2)%nn ! open node B, m2 = # of neighbors

        do k=1,m2 ! loop over neighbor of neighbor 

          if(output(id2)%dist(k).eq.ind2) then ! if 2nd bin satisfied 
            id3=output(id2)%id(k)
            m3=output(id3)%nn ! open node C, m3 = # of neighbors

            do l=1,m3 ! loop of neighbor of neighbor of neighbor 

              if(output(id3)%dist(l).eq.ind3) then
                id4=output(id3)%id(l)
                m4=output(id4)%nn ! open node D, m4 = # of neighbors

                do m=1,m4

                  if(i==output(id4)%id(m) .and. output(id4)%dist(m).eq.ind4) then

                  count=gal(i)+gal(id2)+gal(id3)+gal(id4)

                  if (count==4) then 
                     Zddd(ind1,ind2,ind3,mu1,mu2)=Zddd(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
                     &wgt1(id2)*wgt1(id3)*wgt1(id4))
                  elseif (count==3) then 
                     Zddr(ind1,ind2,ind3,mu1,mu2)=Zddr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
                     &wgt1(id2)*wgt1(id3)*wgt1(id4))
                  elseif (count==2) then
                     Zddr(ind1,ind2,ind3,mu1,mu2)=Zddr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
                     &wgt1(id2)*wgt1(id3)*wgt1(id4))
                  elseif (count==1) then
                     Zdrr(ind1,ind2,ind3,mu1,mu2)=Zdrr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
                     &wgt1(id2)*wgt1(id3)*wgt1(id4))
                  elseif (count==0) then
                     Zrrr(ind1,ind2,ind3,mu1,mu2)=Zrrr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
                     &wgt1(id2)*wgt1(id3)*wgt1(id4))
                  else
                     print*,"weird triplet count"
                  endif
                  exit
              endif
            enddo
           endif
         enddo
       endif 
     enddo
   endif
 enddo
 enddo
!$OMP END PARALLEL DO

end subroutine query_graph_single4

subroutine query_graph_single3(ind1,ind2,ind3)
integer*1 :: ind1,ind2,ind3
integer :: count, m1,m2,m3,id2,id3,m

mu1=1
mu2=1


!$OMP PARALLEL DO schedule(dynamic)  private(i,j,k,l,m,m1,m2,m3,id2,id3,count) & ! , 
!$OMP& shared(wgt1,gal,output,buffer,mu1,mu2,ind1,ind2,ind3)&
!$OMP& reduction(+:Zddd,Zddr,Zdrr,Zrrr)
do i=1,Ndata,1             ! begin loop over all data (and random) point
  if(buffer(i)==1) cycle ! but skip if in buffer region

    m1=output(i)%nn       ! open node A, m1 = # of neighbors

    do j=1,m1             ! loop over neighbor list

      if(output(i)%dist(j).eq.ind1) then !if 1st bin satisfied
        id2=output(i)%id(j)
        m2=output(id2)%nn ! open node B, m2 = # of neighbors

        do k=1,m2 ! loop over neighbor of neighbor 

          if(output(id2)%dist(k).eq.ind2) then ! if 2nd bin satisfied 
            id3=output(id2)%id(k)
            m3=output(id3)%nn ! open node C, m3 = # of neighbors

            do l=1,m3 ! loop of neighbor of neighbor of neighbor 

              if(i==output(id3)%id(l) .and. output(id3)%dist(l).eq.ind3 ) then

                  count=gal(i)+gal(id2)+gal(id3)

                 if (count==3) then 
                     Zddd(ind1,ind2,ind3,mu1,mu2)=Zddd(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
                     &wgt1(id2)*wgt1(id3))
                  elseif (count==2) then
                     Zddr(ind1,ind2,ind3,mu1,mu2)=Zddr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
                     &wgt1(id2)*wgt1(id3))
                  elseif (count==1) then
                     Zdrr(ind1,ind2,ind3,mu1,mu2)=Zdrr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
                     &wgt1(id2)*wgt1(id3))
                  elseif (count==0) then
                     Zrrr(ind1,ind2,ind3,mu1,mu2)=Zrrr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
                     &wgt1(id2)*wgt1(id3))
                  else
                     print*,"weird triplet count"
                  endif
                  exit
           endif
         enddo
       endif 
     enddo
   endif
 enddo
 enddo
!$OMP END PARALLEL DO

end subroutine query_graph_single3

subroutine query_graph_sorted ()
integer nn3,idum,jdum, ninter,count
if(myid==0) print*,'begin querying the graph'

ninter=0

!$OMP PARALLEL DO schedule(dynamic)  private(i,j,nn2,ind1,ind2,ind3,mu1,mu2,nodeid,nn3,idum,jdum,count) & ! , 
!$OMP& shared(wgt1,gal,output,Nrand,myid,buffer)&
!$OMP& reduction(+:Zddd,Zddr,Zdrr,Zrrr,Zdd,Zdr,Zrr)

do i=1,Ndata,1             ! begin loop over all data (and random) point
    if(buffer(i)==1) cycle ! but skip if in buffer region

    nn2=output(i)%nn       ! open node nn2 = # of neighbors
    
    do j=1,nn2             ! loop over neighbor list

        ind1=output(i)%dist(j)
        mu1=output(i)%mu(j)  
        
        count=gal(output(i)%id(j))+gal(i)
        if (count==2) then 
            Zdd(ind1,mu1)=Zdd(ind1,mu1)+wgt1(i)*wgt1(output(i)%id(j))
        elseif (count==1) then 
            Zdr(ind1,mu1)=Zdr(ind1,mu1)+wgt1(i)*wgt1(output(i)%id(j))
        elseif (count==0) then 
            Zrr(ind1,mu1)=Zrr(ind1,mu1)+wgt1(i)*wgt1(output(i)%id(j))
        endif

        nodeid=output(i)%id(j)
        nn3=output(nodeid)%nn

        idum=1
        jdum=1
33      if (idum>nn2 .or. jdum>nn3) then
            goto 34            
        elseif (output(i)%id(idum) < output(nodeid)%id(jdum)) then        
            idum=idum+1
            goto 33
        elseif (output(i)%id(idum) > output(nodeid)%id(jdum)) then
            jdum=jdum+1
            goto 33
        elseif (output(i)%id(idum) == output(nodeid)%id(jdum)) then

          ind2=output(nodeid)%dist(jdum)
          mu2=output(nodeid)%mu(jdum)
          ind3=output(i)%dist(idum) 


!      count=gal(output(i)%id(idum))+gal(nodeid)+gal(i)
      ! if (count==3) then 
      !     Zddd(ind1,ind2,ind3,mu1,mu2)=Zddd(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
      !     &wgt1(output(i)%id(idum))*wgt1(nodeid))
      ! elseif (count==2) then
      !     Zddr(ind1,ind2,ind3,mu1,mu2)=Zddr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
      !     &wgt1(output(i)%id(idum))*wgt1(nodeid))
      ! elseif (count==1) then
      !     Zdrr(ind1,ind2,ind3,mu1,mu2)=Zdrr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
      !     &wgt1(output(i)%id(idum))*wgt1(nodeid))
      ! elseif (count==0) then
      !     Zrrr(ind1,ind2,ind3,mu1,mu2)=Zrrr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
      !     &wgt1(output(i)%id(idum))*wgt1(nodeid))
      ! else
      !   print*,"weird triplet count"
      ! endif


	 if (gal(output(i)%id(idum))==1) then
	   if (gal(nodeid)==1) then
		   if (gal(i)==1) then
       Zddd(ind1,ind2,ind3,mu1,mu2)=Zddd(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
       &wgt1(output(i)%id(idum))*wgt1(nodeid))
     		else
      Zddr(ind1,ind2,ind3,mu1,mu2)=Zddr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
          &wgt1(output(i)%id(idum))*wgt1(nodeid)) 
		  endif
		else
		   if (gal(i)==1) then
      Zddr(ind1,ind2,ind3,mu1,mu2)=Zddr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&  
      &wgt1(output(i)%id(idum))*wgt1(nodeid))
                   else    
      Zdrr(ind1,ind2,ind3,mu1,mu2)=Zdrr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
         &wgt1(output(i)%id(idum))*wgt1(nodeid))
		  endif

		endif

	    else

     if (gal(nodeid)==1) then
       if (gal(i)==1) then
      Zddr(ind1,ind2,ind3,mu1,mu2)=Zddr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
           &wgt1(output(i)%id(idum))*wgt1(nodeid))
        else
      Zdrr(ind1,ind2,ind3,mu1,mu2)=Zdrr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
          &wgt1(output(i)%id(idum))*wgt1(nodeid)) 
      endif
    else
       if (gal(i)==1) then
      Zdrr(ind1,ind2,ind3,mu1,mu2)=Zdrr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
           &wgt1(output(i)%id(idum))*wgt1(nodeid))
                   else    
      Zrrr(ind1,ind2,ind3,mu1,mu2)=Zrrr(ind1,ind2,ind3,mu1,mu2)+(wgt1(i)*&
           &wgt1(output(i)%id(idum))*wgt1(nodeid))
      endif     

	    endif
    endif

            jdum=jdum+1
            idum=idum+1
!            ninter=ninter+1
            goto 33
        endif 
34      continue

    enddo
enddo
!$OMP END PARALLEL DO
!if(myid==0)  print*, '# of triplets ', ninter
deallocate(output)

end subroutine query_graph_sorted

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parseOptions ()
     USE extension, ONLY : getArgument, nArguments
      implicit none

      logical lexist
      integer(kind=4)    :: i,n
      character(len=6)   :: opt
      character(len=2000) :: arg

      n = nArguments()

      if (n < 6 .and. myid==master) then
        print*,' '
        print*,'Not enough input parameters. Please read the following help info'
        print*,' '
        call print_help
        stop
      end if

      i=1

      do while (i<=n) 
         call getArgument(i,opt)

         select case (opt)
         case ('-gal')
            call getArgument(i+1,arg)
            file1 = trim(arg)
            i=i+2
         case ('-ran')
            call getArgument(i+1,arg)
            file2 = trim(arg)
            i=i+2
         case ('-out')
            call getArgument(i+1,arg)
            outfile = trim(arg)
            i=i+2
         case ('-rmin')
            call getArgument(i+1,arg)
            read (arg,*) rmin
            i=i+2
         case ('-rmax')
            call getArgument(i+1,arg)
            read (arg,*) rmax
            i=i+2
         case ('-nbins')
            call getArgument(i+1,arg)
            read (arg,*) nbins
            i=i+2
         case ('-nmu')
            call getArgument(i+1,arg)
            read (arg,*) nmu
            i=i+2
         case ('-rsd')
            RSD=.true.
            i=i+1
         case ('-RSD')
            RSD=.true.
            i=i+1
         case ('-wgt')
              wgt=.true.
              i=i+1
          case ('-cut')
            call getArgument(i+1,arg)
            cut=.true.
            file1 = trim(arg)
            if(myid==master) print*,'Treating input catalogue as subsample'
            if(myid==master) print*,'will ignore -gal -ran options'
            outfile = file1
            i=i+2

         case ('-mpi')
            if(myid==master) print*,'I will read XXX.loadnodes files'
            if(myid==master) print*,'where XXX is the MPI thread number'
            if(myid==master) print*,'make sure # of MPI threads == # loadnode files!!!'
            DOMPI=.true.
            cut=.true.
            i=i+1

         case ('-help')
            call print_help
            stop

         case default
            print '("unknown option ",a6," ignored")', opt
            stop
         end select

      enddo


   end subroutine parseOptions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine print_help
print*,'PURPOSE: Code for calculating the 3PCF of a 3D point set. '
print*,'     '
print*,'     '
print*,'CALLING SEQUENCE:'
print*,'      gramsci [-gal data_file] [-ran random_file] [-out out_file] '
print*,'              [-rmin Rmin] [-rmax Rmax] [-nbins Nbins] '
print*,'              [-cut xxx.loadnodes] [-wgt] [-RSD] [-nmu Nmu] [-mpi]'
print*,' '
print*,'      eg: graph -rmin 10.0 -rmax 12.0 -nmu 10 -nbins 10 -wgt -rsd '
print*,' '
print*,'INPUTS:'
print*,'   '
print*,'       data_file = input catalogue 3-4 columns [x y z (weight)]'
print*,'   '
print*,'       random_file = input random catalogue 3-4 columns [x y z (weight)] - matching above file'
print*,'   '
print*,'       out_file = filename for output correlation func. '
print*,'   '
print*,'       Rmin = Min. radial seperation'
print*,'   '
print*,'       Rmax = Max. radial seperation'
print*,' '
print*,'OPTIONAL:'
print*,'       Nbins = Number of radial bins to calculate [default = 1]'
print*,' '
print*,'       -wgt = if this flag is set, it will use data and random weights '
print*,'              The input gal/ran files should be 4 columns [x y z (weight)]'
print*,' '
print*,'       -RSD = will do anisotropic analysis for redshift distortions and alcock-paczynski effect'
print*,'              In this case the number of angular bins Nmu should be set with -nmu flag'
print*,' '
print*,'       Nmu = Number of angular (mu) bins - linear in range -1<mu<+1 [default = 1]'
print*,' '
print*,'       -cut = Mean the code will work on domain decomposed region files.'
print*,'              These files should have been created with the domain_decompose.sh script'
print*,' '
print*,'       -mpi = I will read XXX.loadnodes with N MPI threads. There should be N loadnode files!'
print*,''
end subroutine print_help
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine count_files ()
implicit none
  Ndata=0
  if(DOMPI) file1=trim(str(myid+1))//'.loadnodes'

  print*, 'opening: ',trim(file1)
  open(11,file=file1,status='unknown')
12 continue
  read(11,*,err=12,end=13)aux
  Ndata=Ndata+1
  goto 12
13 close(11)
  if(myid==0) print*,'Preparing to read ',Ndata, 'data points'
end subroutine count_files

subroutine count_files_2 ()
implicit none
  Ndata=0
  open(11,file=file1,status='unknown')
16 continue
  read(11,*,err=16,end=17)aux
  Ndata=Ndata+1
  goto 16
17 close(11)

  open(11,file=file2,status='unknown')
18 continue
  read(11,*,err=18,end=19)aux
  Ndata=Ndata+1
  goto 18
19 close(11)
  if(myid==0) print*,'Preparing to read ',Ndata, 'data points'
end subroutine count_files_2

subroutine read_files ()
implicit none

if(myid==0)  print*, 'opening ',trim(file1)
open(11,file=file1,status='unknown')
do i=1,Ndata
  if ( wgt ) then
    if ( cut ) then
     read(11,*,end=14)my_array(1:3,i),wgt1(i),gal(i),buffer(i)
    else
     read(11,*,end=14)my_array(1:3,i),wgt1(i),gal(i)
     buffer(i)=0.0
    endif
  else
    if ( cut ) then
      read(11,*,end=14)my_array(1:3,i),gal(i),buffer(i)
      wgt1(i)=1.0
    else
      read(11,*,end=14)my_array(1:3,i),gal(i)
      wgt1(i)=1.0
      buffer(i)=0.0
    endif
 endif
enddo
14 close(11)
if(myid==0)  print*,'Finished reading data file'  
if(myid==0)  print*, 'there are ',sum(gal),' galaxies'
if(myid==0)  print*, 'there are ',Ndata-sum(gal),' randoms'

if(myid==0)  print*, 'there are ',Ndata-sum(buffer),' data points inside buffer'

end subroutine read_files

subroutine read_files_2 ()
implicit none

i=1
if(myid==0)  print*, 'opening ',trim(file1)
open(11,file=file1,status='unknown')
20 continue
  if ( wgt ) then
    read(11,*,err=20,end=21)my_array(1:3,i),wgt1(i)
    gal(i)=1
    buffer(i)=0
  else
    read(11,*,err=20,end=21)my_array(1:3,i)
    gal(i)=1
    wgt1(i)=1.0
    buffer(i)=0
    endif
  i=i+1
goto 20
21 close(11)

if(myid==0)  print*, 'opening ',trim(file2)
open(11,file=file2,status='unknown')
22 continue
  if ( wgt ) then
    if(i>Ndata) goto 23
    read(11,*,err=22,end=23)my_array(1:3,i),wgt1(i)
    gal(i)=0
    buffer(i)=0
  else
    if(i>Ndata) goto 23
    read(11,*,err=22,end=23)my_array(1:3,i)
    gal(i)=0
    wgt1(i)=1.0
    buffer(i)=0
  endif
  i=i+1
goto 22
23 close(11)

if(myid==0)  print*,'Finished reading data file'  
if(myid==0)  print*, 'there are ',sum(gal),' galaxies'
if(myid==0)  print*, 'there are ',Ndata-sum(gal),' randoms'

end subroutine read_files_2



subroutine allocate_arrays ()
  implicit none
  allocate(resultsb(Ndata))

  allocate(Zddd(nbins,nbins,nbins,nmu,nmu))
  allocate(Zddr(nbins,nbins,nbins,nmu,nmu))
  allocate(Zdrr(nbins,nbins,nbins,nmu,nmu))
  allocate(Zrrr(nbins,nbins,nbins,nmu,nmu))
  allocate(crap(nbins,nbins,nbins,nmu,nmu))
  
  allocate(Zdd(nbins,nmu))
  allocate(Zdr(nbins,nmu))
  allocate(Zrr(nbins,nmu))
  allocate(crap2(nbins,nmu))
  
  crap=0.0d0
  Zddd=0.0d0
  Zddr=0.0d0
  Zdrr=0.0d0
  Zrrr=0.0d0
  Zdd=0.0d0
  Zdr=0.0d0
  Zrr=0.0d0
  crap2=0.0d0

end subroutine allocate_arrays

subroutine mpi_collect()
#ifdef MPI
  if(myid==master) crap=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zddd, crap, nbins*nbins*nbins*nmu*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zddd=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zddr, crap, nbins*nbins*nbins*nmu*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zddr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdrr, crap, nbins*nbins*nbins*nmu*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdrr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zrrr, crap, nbins*nbins*nbins*nmu*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zrrr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap2=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdd, crap2, nbins*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdd=crap2
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap2=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdr, crap2, nbins*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdr=crap2
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap2=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zrr, crap2, nbins*nmu, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zrr=crap2
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

#endif

end subroutine

subroutine normalise_counts()
!  real :: avgal,avran
  if ( wgt ) then
    ndd=sum(wgt1,mask = gal .EQ. 1)
    nrr=sum(wgt1,mask = gal .EQ. 0)
    avgal=ndd/sum(gal)
    avran=nrr/(Ndata-sum(gal))
    print*,'weighted # of galaxies / average',ndd, avgal
    print*,'weighted # of randoms / average',nrr, avran
    Zddd=Zddd/(ndd*(ndd-avgal)*(ndd-2*avgal))
    Zrrr=Zrrr/(nrr*(nrr-avran)*(nrr-2*avran))
    Zddr=Zddr/(3.d0*ndd*(ndd-avgal)*nrr)
    Zdrr=Zdrr/(3.d0*ndd*nrr*(nrr-avran))
    
    Zdd=Zdd/(ndd*(ndd-avgal))
    Zdr=Zdr/(2*ndd*nrr)
    Zrr=Zrr/(nrr*(nrr-avran))
else
   ndd=float(Ndata)
   nrr=float(Nrand)
   Zddd=Zddd/(ndd*(ndd-1.)*(ndd-2.))
   Zddr=Zddr/(3*ndd*(ndd-1.)*nrr)
   Zdrr=Zdrr/(3*ndd*nrr*(nrr-1.))
   Zrrr=Zrrr/(nrr*(nrr-1.)*(nrr-2.))
   
    Zdd=Zdd/(ndd*(ndd-1.))
    Zdr=Zdr/(2*ndd*nrr)
    Zrr=Zrr/(nrr*(nrr-1.))
endif

end subroutine

SUBROUTINE dist(vv1,vv2,d)
implicit none
real(kdkind):: vv1(3),vv2(3), d

d=sqrt((vv1(1)-vv2(1))**2 + (vv1(2)-vv2(2))**2 + (vv1(3)-vv2(3))**2)

return
end SUBROUTINE dist

subroutine deallocate_arrays()
implicit none

if (allocated(my_array)) deallocate(my_array)
if (allocated(Zddd)) deallocate(Zddd)
if (allocated(Zddr)) deallocate(Zddr)
if (allocated(Zdrr)) deallocate(Zdrr)
if (allocated(Zrrr)) deallocate(Zrrr)
if (allocated(crap)) deallocate(crap)

if (allocated(Zdd)) deallocate(Zdd)
if (allocated(Zdr)) deallocate(Zdr)
if (allocated(Zrr)) deallocate(Zrr)
if (allocated(crap2)) deallocate(crap2)

if (allocated(wgt1)) deallocate(wgt1)
if (allocated(v1)) deallocate(v1)
if (allocated(v2)) deallocate(v2)
if (allocated(v3)) deallocate(v3)
if (allocated(v4)) deallocate(v4)

if (allocated(resultsb)) deallocate(resultsb)

end subroutine deallocate_arrays


subroutine default_params()

  d=3
  cut=.false.
  wgt=.false.
  nbins=0
  rmin=0.0
  rmax=0.0
  outfile='result.txt'
  RSD=.false.
  nmu=1
  mu1=1
  mu2=1
  
  loadran=.false.
  saveran=.false.

  end subroutine default_params

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

end program  gramsci
