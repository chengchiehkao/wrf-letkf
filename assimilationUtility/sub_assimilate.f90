

!include 'sub_LETKF.f90'

subroutine assimilate(background,analysis,ensembleSize,domain,obs,obsListOfEachGrid)

use derivedType
use basicUtility
use systemUtility
use omp_lib

implicit none

integer,intent(in) :: ensembleSize
type(backgroundInfo),intent(in)  :: background(ensembleSize)
type(backgroundInfo),intent(out) :: analysis(ensembleSize)
type(domainInfo),intent(in)      :: domain(ensembleSize)
type(obsParent) :: obs
type(integerVector),pointer :: obsListOfEachGrid(:,:,:)


real(kind=8),allocatable,dimension(:,:) :: xb_mean,xa_mean
real(kind=8),allocatable,dimension(:,:) :: xb_pert,xa_pert
real(kind=8),allocatable,dimension(:,:) :: yo
real(kind=8),allocatable,dimension(:,:) :: yb
real(kind=8),allocatable,dimension(:,:) :: R

integer obsNumForAssimilation

real(kind=8) :: dh,dz

integer iens,iwe,isn,iz,io,iList  ! loop counter
real ct0,ct1
real(kind=8) wt0,wt1
!================================================

do iens = 1,ensembleSize

    allocate( analysis(iens) % mu( domain(iens)%size_westToEast , domain(iens)%size_southToNorth ) )
    allocate( analysis(iens) % t( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( analysis(iens) % qvapor( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( analysis(iens) % u( domain(iens)%size_westToEast_stag , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( analysis(iens) % v( domain(iens)%size_westToEast , domain(iens)%size_southToNorth_stag , domain(iens)%size_bottomToTop ) )
    allocate( analysis(iens) % w( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop_stag ) )

    analysis(iens) % mu(:,:)       = background(iens) % mu(:,:)  !background(iens) % mu(:,:) + 10.d0
    analysis(iens) % t(:,:,:)      = 0.d0  !background(iens) % t(:,:,:) + 10.d0
    analysis(iens) % qvapor(:,:,:) = 0.d0  !background(iens) % qvapor(:,:,:) + 10.d0
    analysis(iens) % u(:,:,:)      = background(iens) % u(:,:,:)  !background(iens) % u(:,:,:) + 10.d0
    analysis(iens) % v(:,:,:)      = background(iens) % v(:,:,:)  !background(iens) % v(:,:,:) + 10.d0
    analysis(iens) % w(:,:,:)      = background(iens) % w(:,:,:)  !background(iens) % w(:,:,:) + 10.d0
enddo


!$omp parallel do default(private) shared(domain,obsListOfEachGrid,ensembleSize,analysis,background,obs) &
!$omp schedule(dynamic,1)
do iz  = 1,domain(1)%size_bottomToTop
wt0 = omp_get_wtime()
do isn = 1,domain(1)%size_southToNorth
do iwe = 1,domain(1)%size_westToEast

    if ( obsListOfEachGrid(iwe,isn,iz)%vectorSize .eq. 0 ) then  ! no observation means won't update.
        forall(iens=1:ensembleSize) analysis(iens)%t(iwe,isn,iz)      = background(iens)%t(iwe,isn,iz) -300.d0 ! re-substract the offset for wrf
        forall(iens=1:ensembleSize) analysis(iens)%qvapor(iwe,isn,iz) = background(iens)%qvapor(iwe,isn,iz)
        cycle
    else

        obsNumForAssimilation = count( obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(:) ) % available )

        if ( obsNumForAssimilation .eq. 0 ) then
            forall(iens=1:ensembleSize) analysis(iens)%t(iwe,isn,iz)      = background(iens)%t(iwe,isn,iz) -300.d0 ! re-substract the offset for wrf
            forall(iens=1:ensembleSize) analysis(iens)%qvapor(iwe,isn,iz) = background(iens)%qvapor(iwe,isn,iz)
            cycle
        endif

        allocate( xb_mean(2,1) , xb_pert(2,ensembleSize) )
        allocate( xa_mean(2,1) , xa_pert(2,ensembleSize) )

        xb_mean(:,:) = 0.d0
        do iens = 1,ensembleSize
            xb_mean(1,1) = xb_mean(1,1) + background(iens)%t(iwe,isn,iz)
            xb_mean(2,1) = xb_mean(2,1) + background(iens)%qvapor(iwe,isn,iz)
        enddo
        xb_mean(:,:) = xb_mean(:,:) / real(ensembleSize,8)

        forall(iens=1:ensembleSize) xb_pert(1,iens) = background(iens)%t(iwe,isn,iz)      - xb_mean(1,1)
        forall(iens=1:ensembleSize) xb_pert(2,iens) = background(iens)%qvapor(iwe,isn,iz) - xb_mean(2,1)

        allocate( yo(obsNumForAssimilation,1) )
        allocate( yb(obsNumForAssimilation,ensembleSize) )
        allocate( R(obsNumForAssimilation,obsNumForAssimilation) )
        R(:,:) = 0.d0

        io = 0
        do iList = 1,obsListOfEachGrid(iwe,isn,iz)%vectorSize
            if ( obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%available  ) then
                io = io+1
                yo(io,1) = obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%var
                yb(io,:) = obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%background(:)
                dh = greatCircleDistance( &
                  (/domain(1)%lon(iwe,isn),domain(1)%lat(iwe,isn)/) , &
                  (/obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%lon,obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%lat/) &
                  )
                dz = dabs( dlog(domain(1)%pressure(iwe,isn,iz)) - dlog(obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%z) )
                R(io,io) = ( errorFactor(1.d6,0.2d0,dh,dz) * obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%error )**2.d0
            endif
        enddo

        call LETKF( xb_mean(:,:) , xb_pert(:,:) , &
                    xa_mean(:,:) , xa_pert(:,:) , &
                    yo(:,:) , yb(:,:) , R(:,:) , &
                    2 , obsNumForAssimilation , ensembleSize )

        forall(iens=1:ensembleSize) analysis(iens)%t(iwe,isn,iz)      = xa_mean(1,1) + xa_pert(1,iens) -300.d0  ! re-substract the offset for wrf
        forall(iens=1:ensembleSize) analysis(iens)%qvapor(iwe,isn,iz) = xa_mean(2,1) + xa_pert(2,iens)

        deallocate(xb_mean,xa_mean,xb_pert,xa_pert,yo,yb,R)
    endif

enddo
enddo
wt1 = omp_get_wtime()
print*,'Processing: iz=',iz,'walltime=',wt1-wt0,'sec'
enddo
!$omp end parallel do

!================================================
return
stop
end subroutine assimilate
