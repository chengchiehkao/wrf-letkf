
subroutine initializeAnalysis(background,analysis,ensembleSize,domain)

use derivedType

implicit none

integer,intent(in) :: ensembleSize
type(backgroundInfo),intent(in)    :: background(ensembleSize)
type(backgroundInfo),intent(inout) :: analysis(ensembleSize)
type(domainInfo),intent(in)        :: domain(ensembleSize)

integer iens  ! loop counter
!================================================

!$omp parallel do default(shared) &
!$omp private(iens) &
!$omp schedule(dynamic,1)
do iens = 1,ensembleSize

    allocate( analysis(iens) % mu( domain(iens)%size_westToEast , domain(iens)%size_southToNorth ) )
    allocate( analysis(iens) % stratifiedMU( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( analysis(iens) % t( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( analysis(iens) % qvapor( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( analysis(iens) % u( domain(iens)%size_westToEast_stag , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( analysis(iens) % v( domain(iens)%size_westToEast , domain(iens)%size_southToNorth_stag , domain(iens)%size_bottomToTop ) )
    allocate( analysis(iens) % w( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop_stag ) )
    allocate( analysis(iens) % ph( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop_stag ) )

    analysis(iens) % mu(:,:)       = background(iens) % mu(:,:)
    analysis(iens) % stratifiedMU(:,:,:) = background(iens) % stratifiedMU(:,:,:)
    analysis(iens) % t(:,:,:)      = background(iens) % t(:,:,:) -300.d0 ! re-substract the offset for wrf
    analysis(iens) % qvapor(:,:,:) = background(iens) % qvapor(:,:,:)
    analysis(iens) % u(:,:,:)      = background(iens) % u(:,:,:)
    analysis(iens) % v(:,:,:)      = background(iens) % v(:,:,:)
    analysis(iens) % w(:,:,:)      = background(iens) % w(:,:,:)
    analysis(iens) % ph(:,:,:)     = background(iens) % ph(:,:,:)

enddo
!$omp end parallel do


!================================================
return
end subroutine initializeAnalysis
