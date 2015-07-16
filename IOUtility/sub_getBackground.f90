
subroutine getBackground(background,ensembleSize,domain)

use derivedType

implicit none
include 'netcdf.inc'

integer,intent(in)                 :: ensembleSize
type(backgroundInfo),intent(inout) :: background(ensembleSize)
type(domainInfo),intent(in)        :: domain(ensembleSize)

character(len=255) backgroundSrc
character(len=2)  domainSeiralNumInString

integer status,ncID
integer varID_mu , varID_t , varID_qvapor , varID_u , varID_v , varID_w

integer iens
integer iwe,isn,iz
!================================================

write(*,'(a)',advance='no') 'Reading background:'


do iens = 1 , ensembleSize

    write(domainSeiralNumInString,'(i2.2)') iens
    if ( iens.ne.ensembleSize ) write(*,'(1x,a2)',advance='no' ) domainSeiralNumInString
    if ( iens.eq.ensembleSize ) write(*,'(1x,a2)',advance='yes') domainSeiralNumInString
    backgroundSrc = repeat(' ',len(backgroundSrc))
    backgroundSrc = 'input/wrfinput_nc_'//domainSeiralNumInString;

    status = nf_open( backgroundSrc , nf_noWrite , ncID )
    if ( status .ne. nf_noErr ) then
        print*,nf_strError( status )
        stop
    endif

    allocate( background(iens) % mu( domain(iens)%size_westToEast , domain(iens)%size_southToNorth ) )
    allocate( background(iens) % t( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( background(iens) % qvapor( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( background(iens) % u( domain(iens)%size_westToEast_stag , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop ) )
    allocate( background(iens) % v( domain(iens)%size_westToEast , domain(iens)%size_southToNorth_stag , domain(iens)%size_bottomToTop ) )
    allocate( background(iens) % w( domain(iens)%size_westToEast , domain(iens)%size_southToNorth , domain(iens)%size_bottomToTop_stag ) )

    background(iens) % mu(:,:)       = 0.d0
    background(iens) % t(:,:,:)      = 0.d0
    background(iens) % qvapor(:,:,:) = 0.d0
    background(iens) % u(:,:,:)      = 0.d0
    background(iens) % v(:,:,:)      = 0.d0
    background(iens) % w(:,:,:)      = 0.d0


    status = nf_inq_varID( ncID , 'MU' , varID_mu )
    status = nf_get_vara_double( ncID , varID_mu , (/1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,1/) , background(iens)%mu(:,:) )

    status = nf_inq_varID( ncID , 'T' , varID_t )
    status = nf_get_vara_double( ncID , varID_t , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop,1/) , background(iens)%t(:,:,:) )

    status = nf_inq_varID( ncID , 'QVAPOR' , varID_qvapor )
    status = nf_get_vara_double( ncID , varID_qvapor , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop,1/) , background(iens)%qvapor(:,:,:) )

    status = nf_inq_varID( ncID , 'U' , varID_u )
    status = nf_get_vara_double( ncID , varID_u , (/1,1,1,1/) , (/domain(iens)%size_westToEast_stag,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop,1/) , background(iens)%u(:,:,:) )

    status = nf_inq_varID( ncID , 'V' , varID_v )
    status = nf_get_vara_double( ncID , varID_v , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth_stag,domain(iens)%size_bottomToTop,1/) , background(iens)%v(:,:,:) )

    status = nf_inq_varID( ncID , 'W' , varID_w )
    status = nf_get_vara_double( ncID , varID_w , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop_stag,1/) , background(iens)%w(:,:,:) )


    background(iens)%t(:,:,:) = background(iens)%t(:,:,:) + 300.d0  ! Offset of t in WRF is 300K.


    status = nf_close( ncID )
    if ( status .ne. nf_noErr ) then
        print*,nf_strError( status )
        stop
    endif

enddo


!================================================
return
stop
end subroutine getBackground
