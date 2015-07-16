
subroutine convertBackgroundToTemperature(background,ensembleSize,domain)

use derivedType
use basicUtility

implicit none

integer,intent(in)              :: ensembleSize
type(backgroundInfo),intent(inout) :: background(ensembleSize)
type(domainInfo),intent(in)     :: domain(ensembleSize)

integer iens,iwe,isn,iz  ! loop counter

!real(kind=8),external :: convertThetaAndPToTemperature
!================================================


!$omp parallel do default(private) shared(background,domain,ensembleSize) schedule(dynamic,1)
do iens = 1 , ensembleSize

    allocate( background(iens)%normalT( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    background(iens)%normalT(:,:,:) = 0.d0

    do iz = 1 , domain(iens)%size_bottomToTop
    do isn = 1 , domain(iens)%size_southToNorth
    do iwe = 1 , domain(iens)%size_westToEast
        background(iens)%normalT(iwe,isn,iz) = convertThetaAndPToTemperature( background(iens)%t(iwe,isn,iz) , domain(iens)%pressure(iwe,isn,iz) )
    enddo
    enddo
    enddo
enddo
!$omp end parallel do


!================================================
return
stop
end subroutine convertBackgroundToTemperature

