
subroutine getBackground(background,ensembleSize,domain,domainID)

use derivedType

implicit none
include 'netcdf.inc'

integer,intent(in)                 :: ensembleSize
type(backgroundInfo),intent(inout) :: background(ensembleSize)
type(domainInfo),intent(in)        :: domain(ensembleSize)
integer,intent(in)                 :: domainID

character(len=255) backgroundSrc
character(len=2)  domainSeiralNumInString
character(len=2)  domainIDInString

integer ncStatus,ncID
integer varID_mu , varID_t , varID_qvapor , varID_u , varID_v , varID_w , varID_ph , varID_znw , varID_u10 , varID_v10

real(kind=8),allocatable,dimension(:) :: temp_znw , temp_diffZNW

integer iens , iz  ! loop counter
!================================================

write(*,'("D",i2.2,": ",a)',advance='no') domainID,'Reading background:'


do iens = 1 , ensembleSize

    write(domainSeiralNumInString,'(i2.2)') iens
    if ( iens.ne.ensembleSize ) write(*,'(1x,a2)',advance='no' ) domainSeiralNumInString
    if ( iens.eq.ensembleSize ) write(*,'(1x,a2)',advance='yes') domainSeiralNumInString

    write(domainIDInString,'(i2.2)') domainID
    backgroundSrc = repeat(' ',len(backgroundSrc))
    backgroundSrc = 'input/background_d'//domainIDInString//'_'//domainSeiralNumInString;

    ncStatus = nf_open( backgroundSrc , nf_noWrite , ncID )
    if ( ncStatus .ne. nf_noErr ) then
        print*,nf_strError( ncStatus )
        stop
    endif

    if ( iens.ne.ensembleSize ) write(*,'(",")',advance='no' )

    allocate( background(iens) % mu( domain(iens)%size_westToEast , domain(iens)%size_southToNorth ) )
    allocate( background(iens) % u10( domain(iens)%size_westToEast , domain(iens)%size_southToNorth ) )
    allocate( background(iens) % v10( domain(iens)%size_westToEast , domain(iens)%size_southToNorth ) )
    allocate( background(iens) % t( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( background(iens) % qvapor( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( background(iens) % u( domain(iens)%size_westToEast_stag , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( background(iens) % v( domain(iens)%size_westToEast , domain(iens)%size_southToNorth_stag , domain(iens)%size_bottomToTop ) )
    allocate( background(iens) % w( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop_stag ) )
    allocate( background(iens) % ph( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop_stag ) )
    allocate( background(iens) % stratifiedMU( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( temp_znw( domain(iens)%size_bottomToTop_stag ) )
    allocate( temp_diffZNW( domain(iens)%size_bottomToTop ) )

    background(iens) % mu(:,:)       = 0.d0
    background(iens) % u10(:,:)      = 0.d0
    background(iens) % v10(:,:)      = 0.d0
    background(iens) % t(:,:,:)      = 0.d0
    background(iens) % qvapor(:,:,:) = 0.d0
    background(iens) % u(:,:,:)      = 0.d0
    background(iens) % v(:,:,:)      = 0.d0
    background(iens) % w(:,:,:)      = 0.d0
    background(iens) % ph(:,:,:)     = 0.d0
    background(iens) % stratifiedMU(:,:,:) = 0.d0
    temp_znw(:)     = 0.d0
    temp_diffZNW(:) = 0.d0


    ncStatus = nf_inq_varID( ncID , 'MU' , varID_mu )
    ncStatus = nf_get_vara_double( ncID , varID_mu , (/1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,1/) , background(iens)%mu(:,:) )

    ncStatus = nf_inq_varID( ncID , 'U10' , varID_u10 )
    ncStatus = nf_get_vara_double( ncID , varID_u10 , (/1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,1/) , background(iens)%u10(:,:) )

    ncStatus = nf_inq_varID( ncID , 'V10' , varID_v10 )
    ncStatus = nf_get_vara_double( ncID , varID_v10 , (/1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,1/) , background(iens)%v10(:,:) )

    ncStatus = nf_inq_varID( ncID , 'T' , varID_t )
    ncStatus = nf_get_vara_double( ncID , varID_t , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop,1/) , background(iens)%t(:,:,:) )

    ncStatus = nf_inq_varID( ncID , 'QVAPOR' , varID_qvapor )
    ncStatus = nf_get_vara_double( ncID , varID_qvapor , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop,1/) , background(iens)%qvapor(:,:,:) )

    ncStatus = nf_inq_varID( ncID , 'U' , varID_u )
    ncStatus = nf_get_vara_double( ncID , varID_u , (/1,1,1,1/) , (/domain(iens)%size_westToEast_stag,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop,1/) , background(iens)%u(:,:,:) )

    ncStatus = nf_inq_varID( ncID , 'V' , varID_v )
    ncStatus = nf_get_vara_double( ncID , varID_v , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth_stag,domain(iens)%size_bottomToTop,1/) , background(iens)%v(:,:,:) )

    ncStatus = nf_inq_varID( ncID , 'W' , varID_w )
    ncStatus = nf_get_vara_double( ncID , varID_w , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop_stag,1/) , background(iens)%w(:,:,:) )

    ncStatus = nf_inq_varID( ncID , 'PH' , varID_ph )
    ncStatus = nf_get_vara_double( ncID , varID_ph , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop_stag,1/) , background(iens)%ph(:,:,:) )

    ncStatus = nf_inq_varID( ncID , 'ZNW' , varID_znw )
    ncStatus = nf_get_vara_double( ncID , varID_znw , (/1,1/) , (/domain(iens)%size_bottomToTop_stag,1/) , temp_znw(:) )


    temp_diffZNW(:) = temp_znw(1:domain(iens)%size_bottomToTop_stag-1) - temp_znw(2:domain(iens)%size_bottomToTop_stag)

    do iz = 1,domain(iens)%size_bottomToTop
        background(iens) % stratifiedMU(:,:,iz) = temp_diffZNW(iz) * background(iens)%mu(:,:)
    enddo

    background(iens)%t(:,:,:) = background(iens)%t(:,:,:) + 300.d0  ! Offset of t in WRF is 300K.

    ncStatus = nf_close( ncID )
    if ( ncStatus .ne. nf_noErr ) then
        print*,nf_strError( ncStatus )
        stop
    endif

    deallocate(temp_znw,temp_diffZNW)

enddo


!================================================
return
stop
end subroutine getBackground
