

!include 'sub_LETKF.f90'

subroutine assimilate_vGrid(background,analysis,ensembleSize,domain,domain_mean,obs,obsListOfEachGrid,systemParameters)

use derivedType
use basicUtility
use systemUtility
use omp_lib

implicit none

integer,intent(in) :: ensembleSize
type(backgroundInfo),intent(in)  :: background(ensembleSize)
type(backgroundInfo),intent(inout) :: analysis(ensembleSize)
type(domainInfo),intent(in)      :: domain(ensembleSize)
type(domainInfo),intent(in)      :: domain_mean
type(obsParent) :: obs
type(integerVector),pointer :: obsListOfEachGrid(:,:,:)
type(systemParameter),intent(in) :: systemParameters


real(kind=8),allocatable,dimension(:,:) :: xb_mean,xa_mean
real(kind=8),allocatable,dimension(:,:) :: xb_pert,xa_pert
real(kind=8),allocatable,dimension(:,:) :: yo
real(kind=8),allocatable,dimension(:,:) :: yb
real(kind=8),allocatable,dimension(:)   :: R

integer obsNumForAssimilation
real(kind=8) :: dh,dz
integer,parameter :: relaxationZone = 5

integer iens,iwe,isn,iz,io,iList  ! loop counter
real ct0,ct1
real(kind=8) wt0,wt1
!================================================


!$omp parallel do default(private) shared(domain,domain_mean,obsListOfEachGrid,ensembleSize,analysis,background,obs,systemParameters) &
!$omp schedule(dynamic,1)
do iz  = 1,domain_mean%size_bottomToTop
wt0 = omp_get_wtime()
do isn = 1+relaxationZone,domain_mean%size_southToNorth_stag-relaxationZone
do iwe = 1+relaxationZone,domain_mean%size_westToEast-relaxationZone

    if ( obsListOfEachGrid(iwe,isn,iz)%vectorSize .eq. 0 ) then  ! no observation means won't update.
        forall(iens=1:ensembleSize) analysis(iens)%v(iwe,isn,iz) = background(iens)%v(iwe,isn,iz)
        cycle
    else

        obsNumForAssimilation = count( obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(:) ) % available )

        if ( obsNumForAssimilation .eq. 0 ) then
            forall(iens=1:ensembleSize) analysis(iens)%v(iwe,isn,iz) = background(iens)%v(iwe,isn,iz)
            cycle
        endif

        allocate( xb_mean(1,1) , xb_pert(1,ensembleSize) )
        allocate( xa_mean(1,1) , xa_pert(1,ensembleSize) )

        xb_mean(:,:) = 0.d0
        do iens = 1,ensembleSize
            xb_mean(1,1) = xb_mean(1,1) + background(iens)%v(iwe,isn,iz)
        enddo
        xb_mean(:,:) = xb_mean(:,:) / real(ensembleSize,8)

        forall(iens=1:ensembleSize) xb_pert(1,iens) = background(iens)%v(iwe,isn,iz) - xb_mean(1,1)

        allocate( yo(obsNumForAssimilation,1) )
        allocate( yb(obsNumForAssimilation,ensembleSize) )
        allocate( R(obsNumForAssimilation) )
        !R(:,:) = 0.d0

        io = 0
        do iList = 1,obsListOfEachGrid(iwe,isn,iz)%vectorSize
            if ( obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%available ) then
                io = io+1
                yo(io,1) = obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%value
                yb(io,:) = obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%background(:)
                dh = greatCircleDistance( &
                  (/domain_mean%lon_v(iwe,isn),domain_mean%lat_v(iwe,isn)/) , &
                  (/obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%lon,obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%lat/) &
                  )
                dz = dabs( dlog(domain_mean%pressure_v(iwe,isn,iz)) - dlog(obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%z) )
                R(io) = ( errorFactor(systemParameters%rc,systemParameters%rv,dh,dz) * obs%obs( obsListOfEachGrid(iwe,isn,iz)%vector(iList) )%error )**2.d0
            endif
        enddo

        call LETKF( xb_mean(:,:) , xb_pert(:,:) , &
                    xa_mean(:,:) , xa_pert(:,:) , &
                    yo(:,:) , yb(:,:) , R(:) , &
                    1 , obsNumForAssimilation , ensembleSize , systemParameters%inflationFactor )

        forall(iens=1:ensembleSize) analysis(iens)%v(iwe,isn,iz) = xa_mean(1,1) + xa_pert(1,iens)

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
end subroutine assimilate_vGrid
