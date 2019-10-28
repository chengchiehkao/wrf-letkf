
subroutine convertGPHtoGMH(domain,domainSize)

use derivedType
use basicUtility

implicit none

integer,intent(in)             :: domainSize
type(domainInfo),intent(inout) :: domain(domainSize)
real(kind=8),parameter :: invalidValue = -9.d6

integer id,iwe,isn,iz  ! loop conuter
!================================================

do id = 1,domainSize

    allocate( domain(id) % GPH_unstag( domain(id)%size_westToEast , domain(id)%size_southToNorth , domain(id)%size_bottomToTop ) )
    allocate( domain(id) % GMH_unstag( domain(id)%size_westToEast , domain(id)%size_southToNorth , domain(id)%size_bottomToTop ) )

#ifndef PGI
    !$omp parallel do default(private) shared(domain,domainSize,id) schedule(dynamic,1)
#endif
    do isn = 1,domain(id)%size_southToNorth
    do iwe = 1,domain(id)%size_westToEast

        call interp1d( dlog( domain(id) % pressure_w(iwe,isn,domain(id)%size_bottomToTop_stag:1:-1) ) , &
                       domain(id) % GPH(iwe,isn,domain(id)%size_bottomToTop_stag:1:-1) , &
                       domain(id) % size_bottomToTop_stag , &
                       dlog( domain(id) % pressure(iwe,isn,:) ) , &
                       domain(id) % GPH_unstag(iwe,isn,:) , &
                       domain(id) % size_bottomToTop , & 
                       2 , .true. , invalidValue )

    enddo
    enddo
#ifndef PGI
    !$omp end parallel do
#endif

#ifndef PGI
    !$omp parallel do default(private) shared(domain,domainSize,id) schedule(dynamic,1)
#endif
    do iz  = 1,domain(id)%size_bottomToTop
    do isn = 1,domain(id)%size_southToNorth
    do iwe = 1,domain(id)%size_westToEast

        domain(id) % GMH_unstag(iwe,isn,iz) = convertGPHToGMHAtSpecificLatitude( domain(id)%GPH_unstag(iwe,isn,iz) , domain(id)%lat(iwe,isn) )

    enddo
    enddo
    enddo
#ifndef PGI
    !$omp end parallel do
#endif

enddo


!================================================
return
end subroutine convertGPHtoGMH
