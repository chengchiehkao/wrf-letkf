
subroutine getDomain(domain,domainSize)

use derivedType

implicit none
include 'netcdf.inc'
integer,intent(in)           :: domainSize
type(domainInfo),intent(out) :: domain(domainSize)

character(len=255) domainSrc
character(len=2)  domainSeiralNumInString 

integer status,ncID
integer dimID_bottomToTop      , dimID_southToNorth      , dimID_westToEast
integer dimID_bottomToTop_stag , dimID_southToNorth_stag , dimID_westToEast_stag
integer varID_lon , varID_lat , varID_lon_u , varID_lat_u , varID_lon_v , varID_lat_v
integer varID_pressure , varID_GPH , varID_psfc , varID_p_top

real(kind=8),allocatable,dimension(:,:,:) :: pressure_dummy , GPH_dummy
real(kind=8),allocatable,dimension(:,:)   :: psfc_dummy
real(kind=8)                              :: p_top_dummy

integer id  ! loop conuter
!================================================

write(*,'(a)',advance='no') 'Reading domain:'

do id=1,domainSize

    write(domainSeiralNumInString,'(i2.2)') id
    if ( id.ne.domainSize ) write(*,'(1x,a2)',advance='no' ) domainSeiralNumInString
    if ( id.eq.domainSize ) write(*,'(1x,a2)',advance='yes') domainSeiralNumInString
    domainSrc = repeat(' ',len(domainSrc))
    domainSrc = 'input/wrfinput_nc_'//domainSeiralNumInString;

    status = nf_open( domainSrc , nf_noWrite , ncID )
    if ( status .ne. nf_noErr ) then
        print*,nf_strError( status )
        stop
    endif

    status = nf_inq_dimID( ncID , 'bottom_top'       , dimID_bottomToTop       )
    status = nf_inq_dimID( ncID , 'bottom_top_stag'  , dimID_bottomToTop_stag  )
    status = nf_inq_dimID( ncID , 'south_north'      , dimID_southToNorth      )
    status = nf_inq_dimID( ncID , 'south_north_stag' , dimID_southToNorth_stag )
    status = nf_inq_dimID( ncID , 'west_east'        , dimID_westToEast        )
    status = nf_inq_dimID( ncID , 'west_east_stag'   , dimID_westToEast_stag   )

    status = nf_inq_dimlen( ncID , dimID_bottomToTop       , domain(id)%size_bottomToTop       )
    status = nf_inq_dimlen( ncID , dimID_bottomToTop_stag  , domain(id)%size_bottomToTop_stag  )
    status = nf_inq_dimLen( ncID , dimID_southToNorth      , domain(id)%size_southToNorth      )
    status = nf_inq_dimLen( ncID , dimID_southToNorth_stag , domain(id)%size_southToNorth_stag )
    status = nf_inq_dimLen( ncID , dimID_westToEast        , domain(id)%size_westToEast        )
    status = nf_inq_dimLen( ncID , dimID_westToEast_stag   , domain(id)%size_westToEast_stag   )

    allocate( domain(id) % lon( domain(id)%size_westToEast , domain(id)%size_southToNorth ) )
    allocate( domain(id) % lat( domain(id)%size_westToEast , domain(id)%size_southToNorth ) )
    allocate( domain(id) % lon_u( domain(id)%size_westToEast_stag , domain(id)%size_southToNorth ) )
    allocate( domain(id) % lat_u( domain(id)%size_westToEast_stag , domain(id)%size_southToNorth ) )
    allocate( domain(id) % lon_v( domain(id)%size_westToEast , domain(id)%size_southToNorth_stag ) )
    allocate( domain(id) % lat_v( domain(id)%size_westToEast , domain(id)%size_southToNorth_stag ) )
    allocate( domain(id) % pressure ( domain(id)%size_westToEast , domain(id)%size_southToNorth , domain(id)%size_bottomToTop      ) )
    allocate( domain(id) % GPH      ( domain(id)%size_westToEast , domain(id)%size_southToNorth , domain(id)%size_bottomToTop_stag ) )
    allocate( domain(id) % pressure_w ( domain(id)%size_westToEast , domain(id)%size_southToNorth , domain(id)%size_bottomToTop_stag      ) )
    allocate( pressure_dummy ( domain(id)%size_westToEast , domain(id)%size_southToNorth , domain(id)%size_bottomToTop      ) )
    allocate( GPH_dummy      ( domain(id)%size_westToEast , domain(id)%size_southToNorth , domain(id)%size_bottomToTop_stag ) )
    allocate( psfc_dummy     ( domain(id)%size_westToEast , domain(id)%size_southToNorth ) )

    domain(id)%lon(:,:)          = 0.d0
    domain(id)%lat(:,:)          = 0.d0
    domain(id)%lon_u(:,:)        = 0.d0
    domain(id)%lat_u(:,:)        = 0.d0
    domain(id)%lon_v(:,:)        = 0.d0
    domain(id)%lat_v(:,:)        = 0.d0
    domain(id)%pressure(:,:,:)   = 0.d0
    domain(id)%GPH(:,:,:)        = 0.d0
    domain(id)%pressure_w(:,:,:) = 0.d0
    pressure_dummy(:,:,:)        = 0.d0
    GPH_dummy(:,:,:)             = 0.d0

    status = nf_inq_varID( ncID , 'XLONG' , varID_lon )
    status = nf_inq_varID( ncID , 'XLAT'  , varID_lat )
    status = nf_get_vara_double( ncID , varID_lon , (/1,1,1/) , (/domain(id)%size_westToEast,domain(id)%size_southToNorth,1/) , domain(id)%lon(:,:) )
    status = nf_get_vara_double( ncID , varID_lat , (/1,1,1/) , (/domain(id)%size_westToEast,domain(id)%size_southToNorth,1/) , domain(id)%lat(:,:) )

    status = nf_inq_varID( ncID , 'XLONG_U' , varID_lon_u )
    status = nf_inq_varID( ncID , 'XLAT_U'  , varID_lat_u )
    status = nf_get_vara_double( ncID , varID_lon_u , (/1,1,1/) , (/domain(id)%size_westToEast_stag,domain(id)%size_southToNorth,1/) , domain(id)%lon_u(:,:) )
    status = nf_get_vara_double( ncID , varID_lat_u , (/1,1,1/) , (/domain(id)%size_westToEast_stag,domain(id)%size_southToNorth,1/) , domain(id)%lat_u(:,:) )

    status = nf_inq_varID( ncID , 'XLONG_V' , varID_lon_v )
    status = nf_inq_varID( ncID , 'XLAT_V'  , varID_lat_v )
    status = nf_get_vara_double( ncID , varID_lon_v , (/1,1,1/) , (/domain(id)%size_westToEast,domain(id)%size_southToNorth_stag,1/) , domain(id)%lon_v(:,:) )
    status = nf_get_vara_double( ncID , varID_lat_v , (/1,1,1/) , (/domain(id)%size_westToEast,domain(id)%size_southToNorth_stag,1/) , domain(id)%lat_v(:,:) )

    status = nf_inq_varID( ncID , 'P'  , varID_pressure )
    status = nf_get_vara_double( ncID , varID_pressure , (/1,1,1,1/) , (/domain(id)%size_westToEast,domain(id)%size_southToNorth,domain(id)%size_bottomToTop,1/) , pressure_dummy(:,:,:) )
    domain(id)%pressure(:,:,:) = pressure_dummy(:,:,:)
    status = nf_inq_varID( ncID , 'PB' , varID_pressure )
    status = nf_get_vara_double( ncID , varID_pressure , (/1,1,1,1/) , (/domain(id)%size_westToEast,domain(id)%size_southToNorth,domain(id)%size_bottomToTop,1/) , pressure_dummy(:,:,:) )
    domain(id)%pressure(:,:,:) = domain(id)%pressure(:,:,:) + pressure_dummy(:,:,:)

    status = nf_inq_varID( ncID , 'PHB' , varID_GPH )
    status = nf_get_vara_double( ncID , varID_GPH , (/1,1,1,1/) , (/domain(id)%size_westToEast,domain(id)%size_southToNorth,domain(id)%size_bottomToTop_stag,1/) , GPH_dummy(:,:,:) )
    domain(id)%GPH(:,:,:) = GPH_dummy(:,:,:)
    status = nf_inq_varID( ncID , 'PH'  , varID_GPH )
    status = nf_get_vara_double( ncID , varID_GPH , (/1,1,1,1/) , (/domain(id)%size_westToEast,domain(id)%size_southToNorth,domain(id)%size_bottomToTop_stag,1/) , GPH_dummy(:,:,:) )
    domain(id)%GPH(:,:,:) = domain(id)%GPH(:,:,:) + GPH_dummy(:,:,:)
    domain(id)%GPH(:,:,:) = (1.d0/9.81d0) * domain(id)%GPH(:,:,:)

    status = nf_inq_varID( ncID , 'PSFC' , varID_psfc )
    status = nf_get_vara_double( ncID , varID_psfc , (/1,1,1/) , (/domain(id)%size_westToEast,domain(id)%size_southToNorth,1/) , psfc_dummy(:,:) )

    status = nf_inq_varID( ncID , 'P_TOP' , varID_p_top )
    status = nf_get_vara_double( ncID , varID_p_top , (/1/) , (/1/) , p_top_dummy )

    domain(id) % pressure_w(:,:,1) = psfc_dummy(:,:)
    domain(id) % pressure_w(:,:,2:domain(id)%size_bottomToTop_stag-1) = 0.5d0 * ( domain(id)%pressure(:,:,1:domain(id)%size_bottomToTop-1) + domain(id)%pressure(:,:,2:domain(id)%size_bottomToTop) )
    domain(id) % pressure_w(:,:,domain(id)%size_bottomToTop_stag) = p_top_dummy

    status = nf_close( ncID )
    if ( status .ne. nf_noErr ) then
        print*,nf_strError( status )
        stop
    endif

    deallocate( pressure_dummy , GPH_dummy )

enddo

!================================================
return
stop
end subroutine getDomain
