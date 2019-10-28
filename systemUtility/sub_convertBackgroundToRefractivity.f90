
subroutine convertBackgroundToRefractivity(background,ensembleSize,domain)

use derivedType
use basicUtility

implicit none

integer,intent(in)              :: ensembleSize
type(backgroundInfo),intent(inout) :: background(ensembleSize)
type(domainInfo),intent(in)     :: domain(ensembleSize)

integer iens,iwe,isn,iz  ! loop counter

!real(kind=8),external :: convertTAndPAndQvToRefractivity
!================================================


!$omp parallel do default(private) shared(background,domain,ensembleSize) schedule(dynamic,1)
do iens = 1 , ensembleSize

    allocate( background(iens)%refractivity( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    background(iens)%refractivity(:,:,:) = 0.d0

    do iz = 1 , domain(iens)%size_bottomToTop
    do isn = 1 , domain(iens)%size_southToNorth
    do iwe = 1 , domain(iens)%size_westToEast
        background(iens)%refractivity(iwe,isn,iz) = convertTAndPAndQvToRefractivity( background(iens)%normalT(iwe,isn,iz) , domain(iens)%pressure(iwe,isn,iz) , background(iens)%qvapor(iwe,isn,iz) )
    enddo
    enddo
    enddo
enddo
!$omp end parallel do


!================================================
return
end subroutine convertBackgroundToRefractivity
