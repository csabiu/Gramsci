program Ngramsci
  use kdtree2_module
  use kdtree2_precision_module
  use node_module
  use sorting_module
  use iso_fortran_env, only: int8
  implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    INCLUDE "omp_lib.h"
    
  type(node), dimension(:), allocatable :: output
  
  real(kdkind), dimension(:,:), allocatable :: my_array
  real(kdkind), allocatable ::  v1(:),v2(:),v3(:),v4(:),rbin(:)

  type(kdtree2), pointer :: tree
  real(kdkind), allocatable :: N3(:,:,:)
  real(kdkind), allocatable :: N2(:,:,:),wgt1(:),disc(:,:),discN(:,:,:,:),discR(:,:,:,:)
  integer, allocatable ::  buffer(:), bintable(:,:,:,:)
  real :: avgal,avran
  real(kdkind) :: ndd,nrr
  integer :: k, i, j, l,m,n,d,chunk,nmu,nbins,iloop,nrel,cbins
  integer(int8) :: mu1,mu2,mu3,ind1,ind2,ind3,ind4
  character :: file1*2000,file2*2000, ranfile*2000, outfile*2000
  type(kdtree2_result),allocatable :: resultsb(:)
  integer   ::  Ndata,Nrand,nn1,nn2,nnode,nodeid, thread, threads
  real(kdkind)  :: aux,odr,odtheta,theta,start,finish
  real(kdkind) :: rmin, rmax, dr, vec_dist,lrmin
  logical :: wgt, logbins, RSD, loadran,saveran,cut,DOMPI,four_pcf_tetra,four_pcf,three_pcf
  logical :: rancat,four_pcf_box, three_pcf_eq, two_pcf
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

 print*, 'node ',myid, 'of ',ntasks,' checking in...'

!--- set some default parameters ---
  call default_params()

  if(myid==master) print*,'reading options'
  call parseOptions()

  cbins=nbins
  call create_binlookup()
  print*,'Number of binned configurations:',cbins
print*,'number of mu bins:',nmu 

   if(myid==master)  print*,'allocated bins'

   if(cut) then
    call count_files()
  else
    call count_files_2()
    endif


  Allocate(my_array(d,Ndata+Nrand))
  Allocate(wgt1(Ndata+Nrand))
  allocate(buffer(Ndata+Nrand))
  allocate(v1(d))
  allocate(v2(d))
  allocate(v3(d))
  allocate(rbin(nbins+1))
  

print*, "================="
print*, "= defining bins ="
print*, "================="
print*,rmin
  do i=1,nbins+1
    j=i-1
    if(logbins) then
      rbin(i)=10**(log10(rmin)+j*(log10(rmax)-log10(rmin))/float(nbins))
    else
      rbin(i)=rmin+j*(rmax-rmin)/float(nbins)
    endif
    print*, rbin(i)
  enddo 

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
  nn2=kdtree2_r_count_around_point(tp=tree,idxin=j,correltime=-1,r2=rmax*rmax)
  v4(i)=nn2*1.d0
  dr = SUM(v4)/SIZE(v4)
enddo
!print*,maxval(v4),dr,sqrt(SUM((v4-dr)**2)/SIZE(v4))
!dr=dr+2*sqrt(SUM((v4-dr)**2)/SIZE(v4))
!nn2=floor(dr)+1

print*, 'est. memory usage (Gb): ', (Ndata+Nrand)*0.004d0*SUM(v4)/(SIZE(v4)*468911.d0)
print*, 'if this is larger than the RAM in your computer, you should probably kill me'

!!!!!!!!!!!! sampling data !!!!!!!!!!!!!!!
deallocate(v4)
allocate(v4(d))

if(myid==0) print*,'building relationships between nodes'

if(logbins) then
dr=(log10(rmax)-log10(rmin))/nbins
odr=1./dr
lrmin=log10(rmin)
else
dr=(rmax-rmin)/nbins
odr=1./dr
endif

odtheta=(0.5*nmu)


allocate(output(Ndata+Nrand))

!$OMP PARALLEL
!$ threads=OMP_GET_NUM_THREADS()
!$ thread=OMP_GET_THREAD_NUM()
!$OMP END PARALLEL

if(myid==master) print*,'Code running with ',ntasks,' MPI processes each with ',threads,' OMP threads'

call cpu_time(start)
call create_graph(1,999)
call cpu_time(finish)
print '("Creating graph will take ~ ",f10.3," minutes.")',(finish-start)*(Ndata+Nrand)/(60.*1000.*threads)
print*, 'If this takes longer than time to drink a coffee, maybe you should give me more CPUs'
print*, 'Or consider decomposing the domain decomposition option with larger N.'
call create_graph(1000,(Ndata+Nrand))
call cpu_time(finish)
print '("Creating the graph took ",f12.3," seconds.")',(finish-start)/threads


call kdtree2_destroy(tree) 

deallocate(my_array)
deallocate(resultsb) 

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

if(myid==0) print*,'finished building node relationships '


!!!! Timing test !!!!!!
call cpu_time(start)
!if(two_pcf) call query_graph_2pcf(1,3*threads)
!if(three_pcf_eq) call query_graph_equilateral_triangle(1,3*threads)
!!if(four_pcf) call query_4pcf_all(1,3*threads)
!if(four_pcf_tetra) call query_graph_tetrahedron(1,3*threads)
call cpu_time(finish)
print '("Calling graph will take ~ ",f10.3," minutes.")',(finish-start)*(Ndata+Nrand)/(60.*3.*threads*threads)

!!!! Run the desired analysis !!!!!!
call cpu_time(start)
N2=0d0
N3=0d0
if(two_pcf) call query_graph_2pcf(1,Ndata+Nrand)
if(three_pcf) call query_graph_3pcf_all(1,Ndata+Nrand)
if(three_pcf_eq) call query_graph_equilateral_triangle(1,Ndata+Nrand)
if(four_pcf) then
  allocate(discN(cbins,cbins,nbins,6))
  allocate(discR(cbins,cbins,nbins,6))
  discN=0d0
  discR=0d0
  !call query_graph_disc(1,Ndata+Nrand)
  call query_4pcf_all2(1,Ndata+Nrand)
  deallocate(discN)
  deallocate(discR)
endif
if(four_pcf_tetra) then
  allocate(discN(cbins,1,1,6))
  allocate(discR(cbins,1,1,6))
  discN=0d0
  discR=0d0
  call query_graph_tetrahedron(1,Ndata+Nrand)
  deallocate(discN)
  deallocate(discR)
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,size(output)
  call output(i)%destroy()
end do
deallocate(output)

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

!if(myid==master) then
!outfile=trim(outfile)
!open(11,file=outfile,status='unknown')
!do l=1,nbins
!    do k=1,nmu   
!        if(logbins) then
!        write(11,'(8(e14.7,1x))')10**((l-1)*dr+lrmin),10**(l*dr+lrmin),&
!                                 ((float(k)-1.)/odtheta)-1.,(float(k)/odtheta)-1.,&
!             N2(l,k,3),N3(l,k,3),N2(l,k,3)/N3(l,k,3)
!        else
!        write(11,'(8(e14.7,1x))')(l-1)*dr+rmin,l*dr+rmin,((float(k)-1.)/odtheta/2.),(float(k)/odtheta/2.),&
!             N2(l,k,3),N3(l,k,3),N2(l,k,3)/N3(l,k,3)
!        endif
!    enddo
!enddo
!close(11)
!endif

! deallocate the memory for the data.


#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

call deallocate_arrays()

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE( ierr )
#endif

deallocate(buffer)

if(myid==master) then
  print*, "Exit... stage left!"
  stop
else
  stop
endif

contains

subroutine create_graph(istart,iend)
integer :: istart,iend
real :: start, finish

!$OMP PARALLEL DO schedule(dynamic)  private(i,j,k,vec_dist,v1,v2,v3,v4,nn1,nn2,nnode,resultsb,theta) & 
!$OMP& shared(tree,Ndata,Nrand,myid,output,odtheta,nmu,rmin,rmax)
do i=istart,iend

  call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=rmax*rmax,&
       nfound=nn2,nalloc=(Ndata+Nrand),results=resultsb)

  if(rmin>0.0) then
    nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=rmin*rmin)  
  else
    nn1=1
  endif
  nnode=nn2-nn1

  call output(i)%init(nnode)

  v2=my_array(:,i)

  k=0       ! count of valid neighbors added
  do j=nn1+1,nn2
        v1=my_array(:,resultsb(j)%idx)

        call dist(v1,my_array(:,i),vec_dist)

        if(vec_dist<=rmin) cycle
        if(vec_dist>=rmax) cycle
        k=k+1

        if(RSD) then
          v3=0.5*[v2(1)+v1(1),v2(2)+v1(2),v2(3)+v1(3)]
          v4=[v2(1)-v1(1),v2(2)-v1(2),v2(3)-v1(3)] ! connecting vector
          theta=v3(1)*v4(1) + v3(2)*v4(2) + v3(3)*v4(3)
          theta=(theta/((sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)) * (sqrt(v4(1)**2 +v4(2)**2 +v4(3)**2))))

          mu1=floor((theta+1.)*odtheta)+1
          if (mu1 < 1) mu1 = 1
          if (mu1 > nmu) mu1 = nmu
          output(i)%mu(k)=mu1
        else
          output(i)%mu(k)=1
        endif

        if(logbins) then
          output(i)%dist(k)=floor((log10(vec_dist)-lrmin)*odr)+1
        else
          output(i)%dist(k)=floor((vec_dist-rmin)*odr)+1
        endif
        output(i)%id(k)=resultsb(j)%idx
    enddo

    output(i)%nn=k
    if(k>0) call sort2(output(i)%id,output(i)%dist,output(i)%mu,k)

enddo
!$OMP END PARALLEL DO

end subroutine create_graph



subroutine query_graph_2pcf (istart,iend)
integer i,nn2,idum,jdum, ninter,k1,k2,k3,id1,istart,iend
integer(int8) :: ind1,mu
if(myid==0) then
  print*,'begin querying the graph'
  print*,'number of mu bins:',nmu 
endif
!$OMP PARALLEL DO schedule(dynamic)  private(i,k1,nn2,id1,ind1,mu) & ! , 
!$OMP& shared(wgt1,output,Ndata,Nrand,myid,buffer,nmu)&
!$OMP& reduction(+:N2,N3)
do i=istart,iend             ! begin loop over all data points
    if(buffer(i)==1) cycle ! but skip if in buffer region

    nn2=output(i)%nn       ! open node nn2 = # of neighbors

    do k1=1,nn2
      ind1=output(i)%dist(k1)   !edge
      id1= output(i)%id(k1)
      mu=1
      if(RSD) mu=output(i)%mu(k1)

      if(i>Ndata .and. id1>Ndata) then
        N3(ind1,mu,3)=N3(ind1,mu,3)+wgt1(i)*wgt1(id1)
      endif
      N2(ind1,mu,3)=N2(ind1,mu,3)+wgt1(i)*wgt1(id1)

    enddo
enddo
!$OMP END PARALLEL DO

if(myid==master) then
outfile=trim(outfile)
print*,'writing output to: ',  outfile
open(11,file=outfile,status='unknown')
write(11,*) 'r1 min, r1 max, r2 min, r2 max, NN, RR, 2pcf (xi)'

do l=1,nbins
    do k=1,nmu   

        write(11,'(8(e14.7,1x))')rbin(l),rbin(l+1),((float(k)-1.)/odtheta)-1.,(float(k)/odtheta)-1.0,&
             N2(l,k,3),N3(l,k,3),N2(l,k,3)/N3(l,k,3)

    enddo
enddo
close(11)
endif

end subroutine query_graph_2pcf


subroutine query_graph_equilateral_triangle (istart,iend)
integer i,nn2,idum,jdum, ninter,k1,k2,k3,id1,id2,istart,iend
integer(int8) :: ind1,ind2,ind3
if(myid==0) print*,'begin querying the graph'

!$OMP PARALLEL DO schedule(dynamic)  private(i,k1,k2,k3,nn2,id1,id2,ind1,ind2,ind3) & ! , 
!$OMP& shared(wgt1,output,Ndata,Nrand,myid,buffer)&
!$OMP& reduction(+:N2,N3)
do i=istart,iend
!do i=1,Ndata+Nrand,1             ! begin loop over all data (and random) point
    if(buffer(i)==1) cycle ! but skip if in buffer region

    nn2=output(i)%nn       ! open node nn2 = # of neighbors

    do k1=1,nn2
      ind1=output(i)%dist(k1)   !edge
      id1=output(i)%id(k1)

        do k2=k1,nn2
          if(k2==k1)cycle
          ind2=output(i)%dist(k2) !edge
          if(ind1/=ind2) cycle

          id2=output(i)%id(k2)
          call find_dist(id1,id2,ind3) !edge
          if(ind3/=ind1)cycle

          if(RSD) then 
            call find_normal(output(i)%mu(k1),output(i)%mu(k2),ind2)
            N2(ind1,ind2,3)=N2(ind1,ind2,3)+wgt1(i)*wgt1(id1)*wgt1(id2)
          else
            if(i>Ndata .and. id1>Ndata .and. id2>Ndata) then
              N3(ind1,1,3)=N3(ind1,1,3)-wgt1(i)*wgt1(id1)*wgt1(id2)
            endif
            N2(ind1,1,3)=N2(ind1,1,3)+wgt1(i)*wgt1(id1)*wgt1(id2)
          endif

        enddo 
    enddo
enddo
!$OMP END PARALLEL DO


if(myid==master) then
outfile=trim(outfile)
open(11,file=outfile,status='unknown')
do l=1,nbins
    do k=1,nmu   

        write(11,'(8(e14.7,1x))')rbin(l),rbin(l+1),((float(k)-1.)/odtheta/2.),(float(k)/odtheta/2.),&
             N2(l,k,3),N3(l,k,3),N2(l,k,3)/N3(l,k,3)

    enddo
enddo
close(11)
endif

!deallocate(output)

end subroutine query_graph_equilateral_triangle

subroutine query_graph_3pcf_all (istart,iend)
integer i,nn2,idum,jdum, ninter,k1,k2,k3,id1,id2,istart,iend,bin
integer(int8) :: ind1,ind2,ind3

if(myid==0) print*,'Performing 3pcf (all configurations)'
if(myid==0) print*,'begin querying the graph'

!$OMP PARALLEL DO schedule(dynamic)  private(i,k1,k2,k3,nn2,id1,id2,ind1,ind2,ind3,bin) & ! , 
!$OMP& shared(wgt1,output,Ndata,Nrand,myid,buffer)&
!$OMP& reduction(+:N2)&
!$OMP& reduction(+:N3)
do i=istart,iend             ! begin loop over all data (and random) point
    if(buffer(i)==1) cycle ! but skip if in buffer region

    nn2=output(i)%nn       ! open node nn2 = # of neighbors

    do k1=1,nn2
      ind1=output(i)%dist(k1)   !edge
      id1=output(i)%id(k1)

        do k2=k1,nn2
          if(k2==k1)cycle
          ind2=output(i)%dist(k2) !edge
          if(ind2==0) cycle

          id2=output(i)%id(k2)
          call find_dist(id1,id2,ind3) !edge
          if(ind3==0)cycle

          bin=bintable(ind1,ind2,ind3,1)

          if(RSD) then 
            call find_normal(output(i)%mu(k1),output(i)%mu(k2),ind2)
            !N2(bin,ind2,3)=N2(ind1,ind2,3)+wgt1(i)*wgt1(id1)*wgt1(id2)
            if(i>Ndata .and. id1>Ndata .and. id2>Ndata) then
              N3(bin,ind2,3)=N3(bin,ind2,3)-wgt1(i)*wgt1(id1)*wgt1(id2)
            endif
              N2(bin,ind2,3)=N2(bin,ind2,3)+wgt1(i)*wgt1(id1)*wgt1(id2)
          else
            if(i>Ndata .and. id1>Ndata .and. id2>Ndata) then
              N3(bin,1,3)=N3(bin,1,3)-wgt1(i)*wgt1(id1)*wgt1(id2)
            endif
            N2(bin,1,3)=N2(bin,1,3)+wgt1(i)*wgt1(id1)*wgt1(id2)
          endif

        enddo 
    enddo
enddo
!$OMP END PARALLEL DO

if(myid==master) then
outfile=trim(outfile)
open(11,file=outfile,status='unknown')
if(RSD) then
  write(11,*) 'r1 min, r1 max, r2 min, r2 max, r3 min, r3 max, mu min, mu max, NNN, RRR, 3pcf (zeta)'
else
  write(11,*) 'r1 min, r1 max, r2 min, r2 max, r3 min, r3 max, NNN, RRR, 3pcf (zeta)'
endif

cbins=0
do i =1,nbins
  do j=i,nbins
    do k=j,nbins
      cbins=cbins+1
      bintable(i,j,k,1)=cbins

      if(RSD) then
        do l=1,nmu
          write(11,'(11(e14.7,1x))')rbin(i),rbin(i+1),rbin(j),rbin(j+1),rbin(k),rbin(k+1),&
          ((float(l)-1.)/odtheta/2.),(float(l)/odtheta/2.),N2(cbins,l,3),N3(cbins,l,3),N2(cbins,l,3)/N3(cbins,l,3)
        enddo
      else
        write(11,'(9(e14.7,1x))')rbin(i),rbin(i+1),rbin(j),rbin(j+1),rbin(k),rbin(k+1),&
        N2(cbins,1,3),N3(cbins,1,3),N2(cbins,1,3)/N3(cbins,1,3)
      endif

      enddo
    enddo
enddo
close(11)
endif

end subroutine query_graph_3pcf_all


subroutine query_graph_tetrahedron (istart,iend)
integer nn3,idum,jdum, ninter,k1,k2,k3,id1,id2,id3,id4
integer(int8) :: ind1,ind2,ind3,ind4,ind5,ind6,dum
integer istart,iend
if(myid==0) print*,'begin querying the graph'

!$OMP PARALLEL DO schedule(dynamic)  private(dum,i,k1,k2,k3,nn2,id1,id2,id3,id4,ind1,ind2,ind3,ind4,ind5,ind6) & ! , 
!$OMP& shared(wgt1,output,Ndata,Nrand,myid,buffer)&
!$OMP& reduction(+:N2)&
!$OMP& reduction(+:N3)
do i=istart,iend
    if(buffer(i)==1) cycle ! but skip if in buffer region

    nn2=output(i)%nn       ! open node nn2 = # of neighbors

    do k1=1,nn2
      ind1=output(i)%dist(k1)   !edge
      id1=output(i)%id(k1)


        do k2=k1+1,nn2
          ind2=output(i)%dist(k2) !inside
          id2=output(i)%id(k2)
          if(ind1/=ind2) cycle

          call find_dist(id1,id2,ind4) !edge
          if(ind4/=ind1)cycle

            do k3=k2+1,nn2
              ind3=output(i)%dist(k3) !edge
              if(ind3/=ind1)cycle
              id3 =output(i)%id(k3)


              call find_dist(id2,id3,ind5) !edge
              if(ind5/=ind1)cycle

              call find_dist(id1,id3,ind6) !inside
              if(ind1/=ind6) cycle
          
              if(i>Ndata) then
                discN(ind1,1,1,6)=discN(ind1,1,1,6)+wgt1(i)*wgt1(id3)*wgt1(id1)
                if(id3>Ndata .and. id1>Ndata) then
                  discR(ind1,1,1,6)=discR(ind1,1,1,6)+wgt1(i)*wgt1(id3)*wgt1(id1)
                endif
              endif


              if(id3>Ndata) then
                discN(ind1,1,1,1)=discN(ind1,1,1,1)+wgt1(id3)*wgt1(i)*wgt1(id1)
                if(i>Ndata .and. id1>Ndata) then
                  discR(ind1,1,1,1)=discR(ind1,1,1,1)+wgt1(id3)*wgt1(i)*wgt1(id1)
                endif
              endif

              if(i>Ndata) then
                discN(ind1,1,1,4)=discN(ind1,1,1,4)+wgt1(i)*wgt1(id1)*wgt1(id2)
                if(id1>Ndata .and. id2>Ndata) then
                  discR(ind1,1,1,4)=discR(ind1,1,1,4)+wgt1(i)*wgt1(id1)*wgt1(id2)
                endif
              endif

              if(id3>Ndata) then
                discN(ind1,1,1,2)=discN(ind1,1,1,2)+wgt1(id3)*wgt1(i)*wgt1(id2)
                if(i>Ndata .and. id2>Ndata) then
                  discR(ind1,1,1,2)=discR(ind1,1,1,2)+wgt1(id3)*wgt1(i)*wgt1(id2)
                endif
              endif

              if(id1>Ndata) then
                discN(ind1,1,1,3)=discN(ind1,1,1,3)+wgt1(id1)*wgt1(i)*wgt1(id3)
                if(i>Ndata .and. id3>Ndata) then
                  discR(ind1,1,1,3)=discR(ind1,1,1,3)+wgt1(id1)*wgt1(i)*wgt1(id3)
                endif
              endif

              if(i>Ndata) then
                discN(ind1,1,1,5)=discN(ind1,1,1,5)+wgt1(i)*wgt1(id3)*wgt1(id2)
                if(id3>Ndata .and. id2>Ndata) then
                  discR(ind1,1,1,5)=discR(ind1,1,1,5)+wgt1(i)*wgt1(id3)*wgt1(id2)
                endif
              endif

              N2(ind1,1,3)=N2(ind1,1,3)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(ind1,1,3)=N3(ind1,1,3)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif


            enddo
        enddo 
    enddo
enddo
!$OMP END PARALLEL DO

if(myid==master) then
outfile=trim(outfile)
open(11,file=outfile,status='unknown')
do l=1,nbins
    do k=1,nmu

        write(11,'(9(e14.7,1x))')rbin(l),rbin(l+1),((float(k)-1.)/odtheta/2.),(float(k)/odtheta/2.),&
             N2(l,k,3),N3(l,k,3),N2(l,k,3)/N3(l,k,3),&
          (discN(l,1,1,1)/discR(l,1,1,1))*(discN(l,1,1,5)/discR(l,1,1,5))&
        +(discN(l,1,1,3)/discR(l,1,1,3))*(discN(l,1,1,4)/discR(l,1,1,4))&
        +(discN(l,1,1,2)/discR(l,1,1,2))*(discN(l,1,1,6)/discR(l,1,1,6))

    enddo
enddo
close(11)
endif

end subroutine query_graph_tetrahedron

subroutine query_4pcf_all2 (istart,iend)
integer nn3,idum,jdum, ninter,k1,k2,k3,id1,id2,id3,id4
integer(int8) :: ind1,ind2,ind3,ind4,ind5,ind6,dum
integer istart,iend,bin
if(myid==0) print*,'begin querying the graph'

!$OMP PARALLEL DO schedule(dynamic)  private(dum,i,k1,k2,k3,nn2,id1,id2,id3,id4,ind1,ind2,ind3,ind4,ind5,ind6,bin) & ! , 
!$OMP& shared(wgt1,output,Ndata,Nrand,myid,buffer)&
!$OMP& reduction(+:N2,N3,discN,discR)
do i=istart,iend
    if(buffer(i)==1) cycle ! but skip if in buffer region

    nn2=output(i)%nn       ! open node nn2 = # of neighbors
    if(nn2<=3) cycle

    do k1=1,nn2
      ind1=output(i)%dist(k1)   !edge
      id1=output(i)%id(k1)

        do k2=k1+1,nn2
          ind2=output(i)%dist(k2) !inside
          id2=output(i)%id(k2)

          call find_dist(id1,id2,ind4) !edge
          if(ind4==0)cycle

            do k3=k2+1,nn2 
              ind3=output(i)%dist(k3) !edge
              id3 =output(i)%id(k3)

              call find_dist(id2,id3,ind5) !edge
              if(ind5==0) cycle
              call find_dist(id1,id3,ind6) !inside
              if(ind6==0) cycle

              !ind2=1
              !ind6=1

              bin=bintable(ind1,ind4,ind5,ind3)
              if(bin<=0) cycle

              !print*,bin
              N2(bin,ind2,ind6)=N2(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(bin,ind2,ind6)=N3(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif

              !bin=bintable(ind4,ind5,ind3,ind1)
              !N2(bin,ind2,ind6)=N2(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
              !  N3(bin,ind2,ind6)=N3(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !endif

              !bin=bintable(ind5,ind3,ind1,ind4)
              !N2(bin,ind2,ind6)=N2(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
              !  N3(bin,ind2,ind6)=N3(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !endif

              !bin=bintable(ind3,ind1,ind4,ind5)
              !N2(bin,ind2,ind6)=N2(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
              !  N3(bin,ind2,ind6)=N3(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !endif

              !bin=bintable(ind1,ind3,ind5,ind4)
              !N2(bin,ind2,ind6)=N2(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
              !  N3(bin,ind2,ind6)=N3(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !endif

              !bin=bintable(ind3,ind5,ind4,ind1)
              !N2(bin,ind2,ind6)=N2(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
              !  N3(bin,ind2,ind6)=N3(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !endif

              !bin=bintable(ind5,ind4,ind1,ind3)
              !N2(bin,ind2,ind6)=N2(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
              !  N3(bin,ind2,ind6)=N3(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !endif

              !bin=bintable(ind4,ind1,ind3,ind5)
              !N2(bin,ind2,ind6)=N2(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
              !  N3(bin,ind2,ind6)=N3(bin,ind2,ind6)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              !endif
!!!!!!!!!!!!!!!!!
              !goto 1234

              !bin=bintable(ind1,ind4,ind5,ind3)
              if(i>Ndata ) then
                !discN(bin,ind2,ind6,4)=discN(bin,ind2,ind6,4)+abs(wgt1(i))*wgt1(id1)*wgt1(id2)
                !if(id1<Ndata .and. id2<Ndata) then
                  discN(bin,ind2,ind6,4)=discN(bin,ind2,ind6,4)+wgt1(id1)*wgt1(id2)
                !endif
                if(id1>Ndata .and. id2>Ndata) then
                  !discR(bin,ind2,ind6,4)=discR(bin,ind2,ind6,4)+abs(wgt1(i))*wgt1(id1)*wgt1(id2)
                  discR(bin,ind2,ind6,4)=discR(bin,ind2,ind6,4)+wgt1(id1)*wgt1(id2)
                endif
              endif

              if(i>Ndata ) then
                !discN(bin,ind2,ind6,5)=discN(bin,ind2,ind6,5)+abs(wgt1(i))*wgt1(id3)*wgt1(id2)
                !if(id3<Ndata .and. id2<Ndata) then
                  discN(bin,ind2,ind6,5)=discN(bin,ind2,ind6,5)+wgt1(id3)*wgt1(id2)
                !endif
                if(id3>Ndata .and. id2>Ndata) then
                  !discR(bin,ind2,ind6,5)=discR(bin,ind2,ind6,5)+abs(wgt1(i))*wgt1(id3)*wgt1(id2)
                  discR(bin,ind2,ind6,5)=discR(bin,ind2,ind6,5)+wgt1(id3)*wgt1(id2)
                endif
              endif

              if(i>Ndata) then
                !discN(bin,ind2,ind6,6)=discN(bin,ind2,ind6,6)+abs(wgt1(i))*wgt1(id3)*wgt1(id1)
                !if(id3<Ndata .and. id1<Ndata) then
                  discN(bin,ind2,ind6,6)=discN(bin,ind2,ind6,6)+wgt1(id3)*wgt1(id1)
               !endif
                if(id3>Ndata .and. id1>Ndata) then
                  !discR(bin,ind2,ind6,6)=discR(bin,ind2,ind6,6)+abs(wgt1(i))*wgt1(id3)*wgt1(id1)
                  discR(bin,ind2,ind6,6)=discR(bin,ind2,ind6,6)+wgt1(id3)*wgt1(id1)
                endif
              endif

              if(id3>Ndata) then
                !discN(bin,ind2,ind6,2)=discN(bin,ind2,ind6,2)+abs(wgt1(id3))*wgt1(i)*wgt1(id2)
                !if(i<Ndata .and. id2<Ndata) then
                  discN(bin,ind2,ind6,2)=discN(bin,ind2,ind6,2)+wgt1(i)*wgt1(id2)
                !endif
                if(i>Ndata .and. id2>Ndata) then
                  !discR(bin,ind2,ind6,2)=discR(bin,ind2,ind6,2)+abs(wgt1(id3))*wgt1(i)*wgt1(id2)
                  discR(bin,ind2,ind6,2)=discR(bin,ind2,ind6,2)+wgt1(i)*wgt1(id2)
                endif
              endif

              if(id3>Ndata) then
                !discN(bin,ind2,ind6,1)=discN(bin,ind2,ind6,1)+abs(wgt1(id3))*wgt1(i)*wgt1(id1)
                !if(i<Ndata .and. id1<Ndata) then
                  discN(bin,ind2,ind6,1)=discN(bin,ind2,ind6,1)+wgt1(i)*wgt1(id1)
                !endif
                if(i>Ndata .and. id1>Ndata) then
                  !discR(bin,ind2,ind6,1)=discR(bin,ind2,ind6,1)+abs(wgt1(id3))*wgt1(i)*wgt1(id1)
                  discR(bin,ind2,ind6,1)=discR(bin,ind2,ind6,1)+wgt1(i)*wgt1(id1)
                endif
              endif

              if(id1>Ndata) then
                !discN(bin,ind2,ind6,3)=discN(bin,ind2,ind6,3)+abs(wgt1(id1))*wgt1(i)*wgt1(id3)
                !if(i<Ndata .and. id3<Ndata) then
                  discN(bin,ind2,ind6,3)=discN(bin,ind2,ind6,3)+wgt1(i)*wgt1(id3)
                !endif
                if(i>Ndata .and. id3>Ndata) then
                  !discR(bin,ind2,ind6,3)=discR(bin,ind2,ind6,3)+abs(wgt1(id1))*wgt1(i)*wgt1(id3)
                  discR(bin,ind2,ind6,3)=discR(bin,ind2,ind6,3)+wgt1(i)*wgt1(id3)
                endif
              endif
1234 continue

!!!!!!!!!!!!!!!!!
            enddo
        enddo 
    enddo
enddo

!$OMP END PARALLEL DO

if(myid==master) then
print*,'finished counting'
outfile=trim(outfile)
open(11,file=outfile,status='unknown')
write(11,*) 'r1 min, r1 max, r2 min, r2 max, r3 min, r3 max, r4 min, r4 max, NNNN, RRRR, 4pcf (eta)'
cbins=0
do i=1,nbins
  do j=i,nbins
    do k=j,nbins
      do l=k,nbins
        !cbins=cbins+1
        cbins=bintable(i,j,k,l)
        do k1=1,nbins
          do k2=1,nbins
        !k1=1
        !k2=1
        write(11,'(17(e14.7,1x))') rbin(i),rbin(i+1),rbin(j),rbin(j+1),rbin(k),rbin(k+1),rbin(l),rbin(l+1),&
        !rbin(k1),rbin(k1+1),rbin(k2),rbin(k2+1),&
        N2(cbins,k1,k2),N3(cbins,k1,k2),(N2(cbins,k1,k2)/N3(cbins,k1,k2)),&
        (discN(cbins,k1,k2,1)/discR(cbins,k1,k2,1))*(discN(cbins,k1,k2,5)/discR(cbins,k1,k2,5))&
        ,(discN(cbins,k1,k2,3)/discR(cbins,k1,k2,3))*(discN(cbins,k1,k2,4)/discR(cbins,k1,k2,4))&
        ,(discN(cbins,k1,k2,2)/discR(cbins,k1,k2,2))*(discN(cbins,k1,k2,6)/discR(cbins,k1,k2,6))
          enddo
        enddo
        enddo
      enddo
    enddo
enddo
close(11)
endif


end subroutine query_4pcf_all2

subroutine query_4pcf_all (istart,iend)
integer nn3,idum,jdum, ninter,k1,k2,k3,id1,id2,id3,id4
integer(int8) :: ind1,ind2,ind3,ind4,ind5,ind6,dum
integer istart,iend,bin,bin2
real tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
if(myid==0) print*,'begin querying the graph'

!$OMP PARALLEL DO schedule(dynamic)  private(dum,i,k1,k2,k3,nn2,id1,id2,id3,id4,ind1,ind2,ind3,ind4,ind5,ind6,bin,bin2) & ! , 
!$OMP& shared(wgt1,output,Ndata,Nrand,myid,buffer)&
!$OMP& reduction(+:N2,N3,discN,discR)
do i=istart,iend
    if(buffer(i)==1) cycle ! but skip if in buffer region

    nn2=output(i)%nn       ! open node nn2 = # of neighbors

    do k1=1,nn2
      ind1=output(i)%dist(k1)   !edge
      id1=output(i)%id(k1)

        do k2=k1+1,nn2
          ind2=output(i)%dist(k2) !inside
          id2=output(i)%id(k2)

          call find_dist(id1,id2,ind4) !edge
          if(ind4==0)cycle

            do k3=k2+1,nn2 
              ind3=output(i)%dist(k3) !edge
              id3 =output(i)%id(k3)

              call find_dist(id2,id3,ind5) !edge
              if(ind5==0)cycle
              call find_dist(id1,id3,ind6) !inside
              if(ind6==0) cycle


              bin=bintable(ind1,ind2,ind3,1)
              bin=ind1+(ind2-1)*nbins+(ind3-1)*nbins**2+(ind4-1)*nbins**3+(ind5-1)*nbins**4+(ind6-1)*nbins**5
              !bin2=bintable(ind4,ind5,ind6,1)
              bin2=ind4+(ind5-1)*nbins+(ind6-1)*nbins*nbins
              !print*,bin

              call binner(ind1,ind4,ind5,ind3,ind2,ind6,bin)
              N2(bin,1,1)=N2(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(bin,1,1)=N3(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif

              call binner(ind4,ind5,ind3,ind1,ind6,ind2,bin)
              N2(bin,1,1)=N2(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(bin,1,1)=N3(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif

              call binner(ind5,ind3,ind1,ind4,ind2,ind6,bin)
              N2(bin,1,1)=N2(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(bin,1,1)=N3(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif

              call binner(ind3,ind1,ind4,ind5,ind6,ind2,bin)
              N2(bin,1,1)=N2(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(bin,1,1)=N3(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif

              call binner(ind3,ind5,ind4,ind1,ind2,ind6,bin)
              N2(bin,1,1)=N2(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(bin,1,1)=N3(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif

              call binner(ind5,ind4,ind1,ind3,ind6,ind2,bin)
              N2(bin,1,1)=N2(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(bin,1,1)=N3(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif

              call binner(ind4,ind1,ind3,ind5,ind2,ind6,bin)
              N2(bin,1,1)=N2(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(bin,1,1)=N3(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif

              call binner(ind1,ind3,ind5,ind4,ind6,ind2,bin)
              N2(bin,1,1)=N2(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(bin,1,1)=N3(bin,1,1)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif

!!!!!!!!!!!!!!!!!
              !bin=bintable(ind1,ind2,ind3,1)
              if(i>Ndata) then
                discN(ind4,1,1,1)=discN(ind4,1,1,1)+abs(wgt1(i))*wgt1(id1)*wgt1(id2)
                if(id1>Ndata .and. id2>Ndata) then
                  discR(ind4,1,1,1)=discR(ind4,1,1,1)+abs(wgt1(i))*wgt1(id1)*wgt1(id2)
                endif
              endif
              if(i>Ndata) then
                discN(ind5,1,1,1)=discN(ind5,1,1,1)+abs(wgt1(i))*wgt1(id3)*wgt1(id2)
                if(id3>Ndata .and. id2>Ndata) then
                  discR(ind5,1,1,1)=discR(ind5,1,1,1)+abs(wgt1(i))*wgt1(id3)*wgt1(id2)
                endif
              endif
              if(i>Ndata) then
                discN(ind6,1,1,1)=discN(ind6,1,1,1)+abs(wgt1(i))*wgt1(id3)*wgt1(id1)
                if(id3>Ndata .and. id1>Ndata) then
                  discR(ind6,1,1,1)=discR(ind6,1,1,1)+abs(wgt1(i))*wgt1(id3)*wgt1(id1)
                endif
              endif
              if(id3>Ndata) then
                discN(ind2,1,1,1)=discN(ind2,1,1,1)+abs(wgt1(id3))*wgt1(i)*wgt1(id2)
                if(i>Ndata .and. id2>Ndata) then
                  discR(ind2,1,1,1)=discR(ind2,1,1,1)+abs(wgt1(id3))*wgt1(i)*wgt1(id2)
                endif
              endif
              if(id3>Ndata) then
                discN(ind1,1,1,1)=discN(ind1,1,1,1)+abs(wgt1(id3))*wgt1(i)*wgt1(id1)
                if(i>Ndata .and. id1>Ndata) then
                  discR(ind1,1,1,1)=discR(ind1,1,1,1)+abs(wgt1(id3))*wgt1(i)*wgt1(id1)
                endif
              endif
              if(id1>Ndata) then
                discN(ind3,1,1,1)=discN(ind3,1,1,1)+abs(wgt1(id1))*wgt1(i)*wgt1(id3)
                if(i>Ndata .and. id3>Ndata) then
                  discR(ind3,1,1,1)=discR(ind3,1,1,1)+abs(wgt1(id1))*wgt1(i)*wgt1(id3)
                endif
              endif

!!!!!!!!!!!!!!!!!
            enddo
        enddo 
    enddo
enddo
!$OMP END PARALLEL DO

if(myid==master) then
outfile=trim(outfile)
open(11,file=outfile,status='unknown')
cbins=0
do i=1,nbins   !edge
    do j=i,nbins !edge
      do k=j,nbins !edge
        do l=k,nbins !edge
          tmp1=0.0
          tmp2=0.0
          tmp3=0.0
          tmp4=0.0
          tmp5=0.0
          tmp6=0.0
          tmp7=0.0

          do m=1,nbins
          do n=1,nbins
          cbins=cbins+1
        !bin=bintable(i,j,k,1)
        !bin2=bintable(l,m,n,1)
        !bin2=n+(m-1)*nbins+(l-1)*nbins*nbins
        bin=i+(j-1)*nbins+(k-1)*nbins**2+(l-1)*nbins**3+(m-1)*nbins**4+(n-1)*nbins**5

        k1=1
        k2=1

        !tmp1=tmp1+N2(bin,1,1)
        !tmp2=tmp2+N3(bin,1,1)
        !tmp3=tmp3+discN(j,k1,k2,1)
        !tmp4=tmp4+discR(j,k1,k2,1)
        !tmp5=tmp5+discN(n,k1,k2,1)
        !tmp6=tmp6+discR(n,k1,k2,1)

        if(N2(bin,1,k2)==0) cycle
        if(N3(bin,1,k2)==0) cycle

        tmp3=N2(bin,1,k2)/N3(bin,1,k2)

        tmp4=(discN(i,k1,k2,1)/discR(i,k1,k2,1))*(discN(k,k1,k2,1)/discR(k,k1,k2,1))
        tmp5=(discN(j,k1,k2,1)/discR(j,k1,k2,1))*(discN(l,k1,k2,1)/discR(l,k1,k2,1))
        tmp6=(discN(m,k1,k2,1)/discR(m,k1,k2,1))*(discN(n,k1,k2,1)/discR(n,k1,k2,1))
        !print*,tmp3,tmp4,tmp5,tmp6

        if(tmp3==tmp3 .and. abs(tmp3)>0.0 .and. tmp4==tmp4 .and. abs(tmp4)>0.0 .and. &
          tmp5==tmp5 .and. abs(tmp5)>0.0 .and. tmp6==tmp6 .and. abs(tmp6)>0.0  ) then
          tmp2=tmp2+1.
        tmp1=tmp1+(N2(bin,1,k2)/N3(bin,1,k2))-(discN(i,k1,k2,1)/discR(i,k1,k2,1))*(discN(k,k1,k2,1)/discR(k,k1,k2,1))&
        -(discN(j,k1,k2,1)/discR(j,k1,k2,1))*(discN(l,k1,k2,1)/discR(l,k1,k2,1))&
        -(discN(m,k1,k2,1)/discR(m,k1,k2,1))*(discN(n,k1,k2,1)/discR(n,k1,k2,1))
        tmp7=tmp7+(N2(bin,1,k2)/N3(bin,1,k2))
        endif 
        !write(11,'(16(e14.7,1x))') rbin(i),rbin(i+1),rbin(j),rbin(j+1),rbin(k),rbin(k+1),rbin(l),rbin(l+1),&
        !rbin(m),rbin(m+1),rbin(n),rbin(n+1),&
        !N2(bin,1,k2),N3(bin,1,k2),(N2(bin,1,k2)/N3(bin,1,k2)),&
        !(discN(i,k1,k2,1)/discR(i,k1,k2,1))*(discN(k,k1,k2,1)/discR(k,k1,k2,1))&
        !+(discN(j,k1,k2,1)/discR(j,k1,k2,1))*(discN(l,k1,k2,1)/discR(l,k1,k2,1))&
        !+(discN(m,k1,k2,1)/discR(m,k1,k2,1))*(discN(n,k1,k2,1)/discR(n,k1,k2,1))
          enddo
        enddo
        write(11,'(16(e14.7,1x))') rbin(i),rbin(i+1),rbin(j),rbin(j+1),rbin(k),rbin(k+1),rbin(l),rbin(l+1),&
        rbin(m),rbin(m+1),rbin(n),rbin(n+1),&
        tmp1/tmp2,tmp7/tmp2!#,tmp1/tmp2,&
        !(discN(i,k1,k2,1)/discR(i,k1,k2,1))*(discN(m,k1,k2,1)/discR(m,k1,k2,1))&
        !+(tmp3/tmp4)*(tmp5/tmp6)&
        !+(discN(k,k1,k2,1)/discR(k,k1,k2,1))*(discN(l,k1,k2,1)/discR(l,k1,k2,1))

        enddo
      enddo
    enddo
enddo
close(11)
endif

end subroutine query_4pcf_all

subroutine binner(i,j,k,l,m,n,bin)
  integer(int8) :: i,j,k,l,m,n
  integer*4 :: bin
  bin=i+(j-1)*nbins+(k-1)*nbins**2+(l-1)*nbins**3+(m-1)*nbins**4+(n-1)*nbins**5

end subroutine


subroutine query_graph_bipyramid (istart,iend)
integer nn3,idum,jdum, ninter,k1,k2,k3,id1,id2,id3,id4
integer(int8) :: ind1,ind2,ind3,ind4,ind5,ind6,dum
integer istart,iend
if(myid==0) print*,'begin querying the graph'

!$OMP PARALLEL DO schedule(dynamic)  private(dum,i,k1,k2,k3,nn2,id1,id2,id3,id4,ind1,ind2,ind3,ind4,ind5,ind6) & ! , 
!$OMP& shared(wgt1,output,Ndata,Nrand,myid,buffer)&
!$OMP& reduction(+:N2)&
!$OMP& reduction(+:N3)
do i=istart,iend
!do i=1,Ndata+Nrand,1             ! begin loop over all data (and random) point
    if(buffer(i)==1) cycle ! but skip if in buffer region
    
    !do j=1,nbins
    !  N3(j,1,3)=N3(j,1,3)+sum(wgt1(output(i)%id(pack(output(i)%dist(:),output(i)%dist(:)==j))))**2
    !enddo

    !cycle


    nn2=output(i)%nn       ! open node nn2 = # of neighbors

    do k1=1,nn2
      ind1=output(i)%dist(k1)   !edge
      id1=output(i)%id(k1)

        do k2=k1+1,nn2
          !if(k2==k1)cycle
          ind2=output(i)%dist(k2) !inside
          if(ind1/=ind2) cycle

          id2=output(i)%id(k2)
          call find_dist(id1,id2,ind4) !edge
          if(ind4/=ind1)cycle

            do k3=k2+1,nn2
              !if(k3==k2)cycle
              !if(k3==k1)cycle
              ind3=output(i)%dist(k3) !edge
              if(ind3/=ind1)cycle
              id3 =output(i)%id(k3)

              call find_dist(id2,id3,ind5) !edge
              if(ind5/=ind1)cycle
              call find_dist(id1,id3,ind6) !inside
              if(ind1/=ind6) cycle

              N2(ind1,1,3)=N2(ind1,1,3)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              if(i>Ndata .and. id1>Ndata .and. id2>Ndata .and. id3>Ndata ) then
                N3(ind1,1,3)=N3(ind1,1,3)+wgt1(i)*wgt1(id1)*wgt1(id2)*wgt1(id3)
              endif



            enddo
        enddo 
    enddo
enddo
!$OMP END PARALLEL DO

if(myid==master) then
outfile=trim(outfile)
open(11,file=outfile,status='unknown')
do l=1,nbins
    do k=1,nmu   
        if(logbins) then
        write(11,'(8(e14.7,1x))')10**((l-1)*dr+lrmin),10**(l*dr+lrmin),&
                                 ((float(k)-1.)/odtheta)-1.,(float(k)/odtheta)-1.,&
             N2(l,k,3),N3(l,k,3),N2(l,k,3)/N3(l,k,3)
        else
        write(11,'(8(e14.7,1x))')(l-1)*dr+rmin,l*dr+rmin,((float(k)-1.)/odtheta/2.),(float(k)/odtheta/2.),&
             N2(l,k,3),N3(l,k,3),N2(l,k,3)/N3(l,k,3)
        endif
    enddo
enddo
close(11)
endif

end subroutine query_graph_bipyramid


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
            rancat=.true.
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
            if(nmu>=2) then
              RSD=.true.
              if(myid==master) print*,'Anisotropic analysis requested'
            i=i+2
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
            DOMPI=.true.
            i=i+1

          case ('-log')
            logbins=.true.
            if(myid==master) print*,'Using logarithmically spaced bins'
            i=i+1

         case ('-help')
            call print_help
            stop

          case ('-tetra')
            four_pcf_tetra=.true.
            i=i+1

          case ('-3pcf')
            three_pcf=.true.
            i=i+1

          case ('-4pcf')
            four_pcf=.true.
            i=i+1
          case ('-equi')
            three_pcf_eq=.true.
            i=i+1

          case ('-2pcf')
            two_pcf=.true.
            i=i+1

         case default
            print '("unknown option ",a," ignored")', trim(opt)
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
print*,'      Ngramsci[-gal galaxy file][-ran ranfile (optional)]'
print*,'           [-rmin Rmin] [-rmax Rmax] [-nbins Nbins] '
print*,'           [-RSD RSD] [-nmu Nmu] [-out out_file]'
print*,' '
print*,'      eg: Ngramsci -rmin 10.0 -rmax 12.0 -nbins 10 [choose one: -2pcf -3pcf -4pcf'
print*,' '
print*,'INPUTS:'
print*,'   '
print*,'       out_file = filename for output correlation func. '
print*,'   '
print*,'       Rmin = Min. radial seperation'
print*,'   '
print*,'       Rmax = Max. radial seperation'
print*,' '
print*,'OPTIONAL:'
print*,'       Nbins = Number of radial bins to calculate'
print*,' '
print*,'       Nmu = Number of angular (mu) bins - linear in range -1<mu<+1'
print*,' '
print*,'       RSD = logical value to request anisotropic analysis.'
print*,'              In this case the number of angular bins Nmu should be set'
print*,' '
print*,''
end subroutine print_help
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine count_files ()
implicit none
  Ndata=0
  if(DOMPI) file1=trim(str(myid+1))//'.loadnodes'

  print*, file1
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
  Nrand=0
  open(11,file=file1,status='unknown')
16 continue
  read(11,*,err=16,end=17)aux
  Ndata=Ndata+1
  goto 16
17 close(11)

if(rancat) then
  open(11,file=file2,status='unknown')
18 continue
   read(11,*,err=18,end=19)aux
   Nrand=Nrand+1
   goto 18
 19 close(11)
endif
if(myid==0) print*,'Preparing to read ',Ndata, 'data points'
if(myid==0) print*,'Preparing to read ',Nrand, 'random points'
end subroutine count_files_2

subroutine read_files ()
implicit none

if(myid==0)  print*, 'opening ',trim(file1)
open(11,file=file1,status='unknown')
do i=1,Ndata
    if ( cut ) then
     read(11,*,end=14)my_array(1:3,i),wgt1(i),buffer(i)
    else
     read(11,*,end=14)my_array(1:3,i),wgt1(i)
     buffer(i)=0.0
    endif
enddo
14 close(11)
if(myid==0)  print*,'Finished reading data file'  
!if(myid==0)  print*, 'there are ',sum(gal),' galaxies'
!if(myid==0)  print*, 'there are ',Ndata-sum(gal),' randoms'


if(myid==0)  print*, 'there are ',Ndata-sum(buffer),' data points inside buffer'
end subroutine read_files

subroutine read_files_2 ()
implicit none

if(myid==0)  print*, 'opening ',trim(file1)
open(11,file=file1,status='unknown')
!20 continue
20 do i=1,Ndata
  read(11,*,err=20,end=21)my_array(1:3,i),wgt1(i)
  buffer(i)=0

enddo 
21 close(11)

if(rancat) then
if(myid==0)  print*, 'opening ',trim(file2)
open(11,file=file2,status='unknown')

22 do i=Ndata+1,Ndata+Nrand
  read(11,*,err=22,end=23)my_array(1:3,i),wgt1(i)
  buffer(i)=0
enddo
23 close(11)
endif

!normalise!!

wgt1(1:Ndata)=wgt1(1:Ndata)/sum(wgt1(1:Ndata))
wgt1(Ndata+1:Ndata+Nrand)=-1.d0*wgt1(Ndata+1:Ndata+Nrand)/sum(wgt1(Ndata+1:Ndata+Nrand))

if(myid==0)  print*,'Finished reading data file'  
if(myid==0)  print*, 'sum of weights: ',sum(wgt1)

end subroutine read_files_2

subroutine create_binlookup()
integer i,j,k,l,m,n


if(four_pcf) then
allocate(bintable(nbins,nbins,nbins,nbins))
cbins=0
!do i =1,nbins
!  do j=i,nbins
!    do k=j,nbins
!      do l=k,nbins
do i =1,nbins
  do j=i,nbins
    do k=j,nbins
      do l=k,nbins
        cbins=cbins+1
        bintable(i,j,k,l)=cbins

        !bintable(i,j,l,k)=cbins
        !bintable(i,k,j,l)=cbins
        !bintable(i,k,l,j)=cbins
        !bintable(i,l,j,k)=cbins
        !bintable(i,l,k,j)=cbins
        !bintable(j,i,k,l)=cbins
        !bintable(j,i,l,k)=cbins
        !bintable(j,k,i,l)=cbins
        !bintable(j,k,l,i)=cbins
        !bintable(j,l,i,k)=cbins
        !bintable(j,l,k,i)=cbins
        !bintable(k,i,j,l)=cbins
        !bintable(k,i,l,j)=cbins
        !bintable(k,j,i,l)=cbins
        !bintable(k,j,l,i)=cbins
        !bintable(k,l,i,j)=cbins
        !bintable(k,l,j,i)=cbins
        !bintable(l,i,j,k)=cbins
        !bintable(l,i,k,j)=cbins
        !bintable(l,j,i,k)=cbins
        !bintable(l,j,k,i)=cbins
        !bintable(l,k,i,j)=cbins
        !bintable(l,k,j,i)=cbins
      enddo
    enddo
  enddo
enddo

!cbins=0
!do i =1,nbins
!  do j=i,nbins
!    do k=j,nbins
!      cbins=cbins+1
!      bintable(i,j,k,1)=cbins
!      bintable(i,k,j,1)=cbins
!      bintable(j,i,k,1)=cbins
!      bintable(j,k,i,1)=cbins
!      bintable(k,i,j,1)=cbins
!      bintable(k,j,i,1)=cbins
!    enddo
!  enddo
!enddo

elseif(three_pcf) then
  allocate(bintable(nbins,nbins,nbins,1))
cbins=0
do i =1,nbins
  do j=i,nbins
    do k=j,nbins
      cbins=cbins+1
      bintable(i,j,k,1)=cbins
      bintable(i,k,j,1)=cbins
      bintable(j,i,k,1)=cbins
      bintable(j,k,i,1)=cbins
      bintable(k,i,j,1)=cbins
      bintable(k,j,i,1)=cbins
    enddo
  enddo
enddo
endif

end subroutine create_binlookup

subroutine allocate_arrays ()
  implicit none
  allocate(resultsb(Ndata+Nrand))

  !allocate(N3(nbins,nbins,nbins,nmu,nmu,4))
  !allocate(N3_tmp(nbins,nbins,nbins,nmu,nmu,4))
  
  if(four_pcf) then 
    !allocate(N2(nbins**6,1,1))
    !allocate(N3(nbins**6,1,1))
    allocate(N2(cbins,cbins,cbins))
    allocate(N3(cbins,cbins,cbins))
  else
    allocate(N2(cbins,nmu,3))
    allocate(N3(cbins,nmu,3))
  endif
  !allocate(N2_tmp(cbins,nmu,3))
  
  !N3_tmp=0.0d0
  N3=0.0d0

  N2=0.0d0
  !N2_tmp=0.0d0

end subroutine allocate_arrays

subroutine mpi_collect()
#ifdef MPI
!if(myid==master) N3_tmp=0.d0
!call MPI_Barrier(MPI_COMM_WORLD,ierr)
!call MPI_REDUCE( N3, N3_tmp, nbins*nbins*nbins*nmu*nmu*4, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
!if ( myid == master ) then !in master thread
!N3=N3_tmp
!endif
!call MPI_Barrier(MPI_COMM_WORLD,ierr)

!if(myid==master) N2_tmp=0.d0
!call MPI_Barrier(MPI_COMM_WORLD,ierr)
!call MPI_REDUCE( N2, N2_tmp, nbins*nmu*3, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
!if ( myid == master ) then !in master thread
!N2=N2_tmp
!endif
!call MPI_Barrier(MPI_COMM_WORLD,ierr)

#endif

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
!if (allocated(N3)) deallocate(N3)
!if (allocated(N3_tmp)) deallocate(N3_tmp)

if (allocated(N2)) deallocate(N2)
!if (allocated(N2_tmp)) deallocate(N2_tmp)

if (allocated(wgt1)) deallocate(wgt1)
if (allocated(v1)) deallocate(v1)
if (allocated(v2)) deallocate(v2)
if (allocated(v3)) deallocate(v3)
if (allocated(v4)) deallocate(v4)

if (allocated(resultsb)) deallocate(resultsb)

end subroutine deallocate_arrays

subroutine find_dist2(id1,id2,ind)
implicit none
integer i
INTEGER, INTENT(IN)                :: id1, id2
integer(int8), INTENT(OUT)                :: ind

if(output(id1)%id(1)>id2 .or. output(id1)%id(output(id1)%nn)<id2) then
  ind=0
  return
endif

ind=0
  do i=1,output(id1)%nn
    if(output(id1)%id(i)==id2) then
      ind=output(id1)%dist(i)
      return
    endif
  enddo 
end subroutine



subroutine find_dist(id1,id2, ind) !Binary chopper. Find i such that X = A(i)
  implicit none
  integer i
  INTEGER, INTENT(IN)                :: id1, id2
  integer(int8), INTENT(OUT)                :: ind
  INTEGER L,R,P    !Fingers.

        
  if(id1<=0) then
          print*, 'ERROR: invalid node index', id1
        endif

        if(output(id1)%id(1)>id2 .or. output(id1)%id(output(id1)%nn)<id2) then
          ind=0
          return
        endif

        L = 0     !Establish outer bounds, to search A(L+1:R-1).
        R = output(id1)%nn + 1   !L = first - 1; R = last + 1.
    1   P = (R - L)/2   !Probe point. Beware INTEGER overflow with (L + R)/2.
        IF (P.LE.0) GO TO 5 !Aha! Nowhere!! The span is empty.
        P = P + L   !Convert an offset from L to an array index.
        IF (id2 - output(id1)%id(P)) 3,4,2 !Compare to the probe point.
    2   L = P     !A(P) < X. Shift the left bound up: X follows A(P).
        GO TO 1     !Another chop.
    3   R = P     !X < A(P). Shift the right bound down: X precedes A(P).
        GO TO 1     !Try again.
    4   ind = output(id1)%dist(P)   !A(P) = X. So, X is found, here!
       RETURN     !Done.
    5   ind = 0   !X is not found. Insert it at L + 1, i.e. at A(1 - FINDI).
      END SUBROUTINE  !A's values need not be all different, merely in order. 

subroutine find_normal(mu1,mu2,mun)
implicit none
integer(int8), INTENT(in)                :: mu1,mu2
integer(int8), INTENT(OUT)                :: mun
real :: mu11,mu22

  mu11=((mu1-0.5)/odtheta) -1.0 !+ 1./nmu
  mu22=((mu2-0.5)/odtheta) -1.0 !+ 1./nmu
  
  mun=floor((1.1547*(0.75 - mu11*mu11 - mu22*mu22 + (mu11*mu22))**0.5)*nmu)+1
  if (mun < 1) mun = 1
  if (mun > nmu) mun = nmu
  !print*,mu11,mu22,((.75 - mu11*mu11 - mu22*mu22 + (mu11*mu22))),mun
  !print*,mu1,mu11,mu2,mu22,(0.75 - mu11*mu11 - mu22*mu22 + mu11*mu22),1.1547*(0.75 - mu11*mu11 - mu22*mu22 + mu11*mu22)**0.5,mun
end subroutine find_normal

subroutine default_params()

  d=3
  cut=.false.
  wgt=.true.
  logbins=.false.
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
  rancat=.false.

  four_pcf_tetra=.false.
  four_pcf_box=.false.
  three_pcf_eq=.false.
  three_pcf=.false.
  four_pcf=.false.
  two_pcf=.false.

end subroutine default_params

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

end program  Ngramsci
