
subroutine getMeanDomain(domain,ensembleSize,domain_mean)

use derivedType

implicit none
integer,intent(in)          :: ensembleSize
type(domainInfo),intent(in) :: domain(ensembleSize)
type(domainInfo)            :: domain_mean

integer iens  ! loop counter
!================================================

!
!  Assume dimensions of each domain are all the same.
!
allocate( domain_mean % lon( domain(1)%size_westToEast , domain(1)%size_southToNorth ) )
allocate( domain_mean % lat( domain(1)%size_westToEast , domain(1)%size_southToNorth ) )
allocate( domain_mean % lon_u( domain(1)%size_westToEast_stag , domain(1)%size_southToNorth ) )
allocate( domain_mean % lat_u( domain(1)%size_westToEast_stag , domain(1)%size_southToNorth ) )
allocate( domain_mean % lon_v( domain(1)%size_westToEast , domain(1)%size_southToNorth_stag ) )
allocate( domain_mean % lat_v( domain(1)%size_westToEast , domain(1)%size_southToNorth_stag ) )
allocate( domain_mean % pressure ( domain(1)%size_westToEast , domain(1)%size_southToNorth , domain(1)%size_bottomToTop      ) )
allocate( domain_mean % GPH      ( domain(1)%size_westToEast , domain(1)%size_southToNorth , domain(1)%size_bottomToTop_stag ) )
allocate( domain_mean % pressure_w ( domain(1)%size_westToEast , domain(1)%size_southToNorth , domain(1)%size_bottomToTop_stag      ) )

!
!  Assume horizontal coordinates of each domain are all the same.
!
domain_mean % lon(:,:)   = domain(1) % lon(:,:)
domain_mean % lat(:,:)   = domain(1) % lat(:,:)
domain_mean % lon_u(:,:) = domain(1) % lon_u(:,:)
domain_mean % lat_u(:,:) = domain(1) % lat_u(:,:)
domain_mean % lon_v(:,:) = domain(1) % lon_v(:,:)
domain_mean % lat_v(:,:) = domain(1) % lat_v(:,:)


domain_mean % pressure(:,:,:)   = 0.d0
domain_mean % GPH(:,:,:)        = 0.d0
domain_mean % pressure_w(:,:,:) = 0.d0

do iens = 1,ensembleSize
    domain_mean % pressure(:,:,:)   = domain_mean % pressure(:,:,:)   + domain(iens) % pressure(:,:,:)
    domain_mean % GPH(:,:,:)        = domain_mean % GPH(:,:,:)        + domain(iens) % GPH(:,:,:)
    domain_mean % pressure_w(:,:,:) = domain_mean % pressure_w(:,:,:) + domain(iens) % pressure_w(:,:,:)
enddo

domain_mean % pressure(:,:,:)   = (1d0/real(ensembleSize,8)) * domain_mean % pressure(:,:,:)
domain_mean % GPH(:,:,:)        = (1d0/real(ensembleSize,8)) * domain_mean % GPH(:,:,:)
domain_mean % pressure_w(:,:,:) = (1d0/real(ensembleSize,8)) * domain_mean % pressure_w(:,:,:)


!================================================
return
stop
end subroutine getMeanDomain
