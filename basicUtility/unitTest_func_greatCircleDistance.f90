

implicit none

integer n
real(kind=8),allocatable,dimension(:,:) :: grid,obs,gridInRad,obsInRad
real(kind=8),allocatable,dimension(:) :: sin_gridInRad,sin_obsInRad
real(kind=8),allocatable,dimension(:) :: cos_gridInRad,cos_obsInRad
real(kind=8),allocatable,dimension(:)   :: dist
real(kind=8),external :: greatCircleDistance,greatCircleDistance_preCalc
real(kind=8),parameter :: radian = 57.29577951308232d0
integer i,j

real t0,t1
!================================================

print*,'n='
read(*,*) n

allocate(grid(2,n),obs(2,300),dist(n))
allocate(gridInRad(2,n),obsInRad(2,300))
allocate(sin_gridInRad(n),sin_obsInRad(n))
allocate(cos_gridInRad(n),cos_obsInRad(n))
call random_seed()
call random_number(grid(:,:))
call random_number(obs(:,:))

grid(:,1) = 90.d0 + 50.d0*grid(:,1)
grid(:,2) =  8.d0 + 35.d0*grid(:,2)
obs(:,1) = 90.d0 + 50.d0*obs(:,1)
obs(:,2) =  8.d0 + 35.d0*obs(:,2)


call cpu_time(t0)
do i=1,n
    do j=1,300
        dist(i) = greatCircleDistance(grid(:,i),obs(:,j))
    enddo
enddo
call cpu_time(t1)
print*,'cpu time=',t1-t0,'sec'


call cpu_time(t0)
gridInRad(:,:) = gridInRad(:,:)/radian
obsInRad(:,:)  = obsInRad(:,:) /radian
sin_gridInRad(:) = dsin(sin_gridInRad(:))
sin_obsInRad(:)  = dsin(sin_obsInRad(:))
cos_gridInRad(:) = dcos(cos_gridInRad(:))
cos_obsInRad(:)  = dcos(cos_obsInRad(:))
do i=1,n
    do j=1,300
        dist(i) = greatCircleDistance_preCalc(gridInRad(1,i),obsInRad(1,j),sin_gridInRad(i),sin_obsInRad(j), &
        cos_gridInRad(i),cos_obsInRad(j))
    enddo
enddo
call cpu_time(t1)
print*,'cpu time=',t1-t0,'sec'


!================================================
stop
end


