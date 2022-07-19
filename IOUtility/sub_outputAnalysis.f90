
subroutine outputAnalysis(analysis,ensembleSize,domain,domainID)

use derivedType

implicit none
include 'netcdf.inc'

integer,intent(in)                 :: ensembleSize
type(backgroundInfo),intent(inout) :: analysis(ensembleSize)
type(domainInfo),intent(in)        :: domain(ensembleSize)
integer,intent(in)                 :: domainID

character(len=255) backgroundSrc,analysisSrc
character(len=2)  domainSeiralNumInString
character(len=2)  domainIDInString

integer ncStatus,ncID
integer varID_mu , varID_t , varID_qvapor , varID_u , varID_v , varID_w ,varID_ph

integer iens
integer iwe,isn,iz
!================================================


do iens =1,ensembleSize

    analysis(iens)%mu(:,:) = sum( analysis(iens)%stratifiedMU(:,:,:) , 3 )

    where ( analysis(iens)%qvapor(:,:,:) .lt. 0.d0 )
        analysis(iens)%qvapor = 0.d0
    end where

enddo

write(*,'("D",i2.2,": ",a)',advance='no') domainID,'Output analysis:'


do iens = 1 , ensembleSize

    write(domainSeiralNumInString,'(i2.2)') iens
    if ( iens.ne.ensembleSize ) write(*,'(1x,a2)',advance='no' ) domainSeiralNumInString
    if ( iens.eq.ensembleSize ) write(*,'(1x,a2)',advance='yes') domainSeiralNumInString

    write(domainIDInString,'(i2.2)') domainID
    analysisSrc = repeat(' ',len(analysisSrc))
    analysisSrc = 'output/analysis_d'//domainIDInString//'_'//domainSeiralNumInString;
    backgroundSrc = repeat(' ',len(backgroundSrc))
    backgroundSrc = 'input/background_d'//domainIDInString//'_'//domainSeiralNumInString;

    call system('cp '//backgroundSrc//' '//analysisSrc)


    ncStatus = nf_open( analysisSrc , nf_write , ncID )
    if ( ncStatus .ne. nf_noErr ) then
        print*,nf_strError( ncStatus )
        stop
    endif

    if ( iens.ne.ensembleSize ) write(*,'(",")',advance='no' )

    ncStatus = nf_inq_varID( ncID , 'MU' , varID_mu )
    ncStatus = nf_put_vara_double( ncID , varID_mu , (/1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,1/) , analysis(iens)%mu(:,:) )
 
    ncStatus = nf_inq_varID( ncID , 'T' , varID_t )
    ncStatus = nf_put_vara_double( ncID , varID_t , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop,1/) , analysis(iens)%t(:,:,:) )
 
    ncStatus = nf_inq_varID( ncID , 'QVAPOR' , varID_qvapor )
    ncStatus = nf_put_vara_double( ncID , varID_qvapor , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop,1/) , analysis(iens)%qvapor(:,:,:) )
 
    ncStatus = nf_inq_varID( ncID , 'U' , varID_u )
    ncStatus = nf_put_vara_double( ncID , varID_u , (/1,1,1,1/) , (/domain(iens)%size_westToEast_stag,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop,1/) , analysis(iens)%u(:,:,:) )
 
    ncStatus = nf_inq_varID( ncID , 'V' , varID_v )
    ncStatus = nf_put_vara_double( ncID , varID_v , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth_stag,domain(iens)%size_bottomToTop,1/) , analysis(iens)%v(:,:,:) )
 
    ncStatus = nf_inq_varID( ncID , 'W' , varID_w )
    ncStatus = nf_put_vara_double( ncID , varID_w , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop_stag,1/) , analysis(iens)%w(:,:,:) )
 
    ncStatus = nf_inq_varID( ncID , 'PH' , varID_ph )
    ncStatus = nf_put_vara_double( ncID , varID_ph , (/1,1,1,1/) , (/domain(iens)%size_westToEast,domain(iens)%size_southToNorth,domain(iens)%size_bottomToTop_stag,1/) , analysis(iens)%ph(:,:,:) )

    ncStatus = nf_close( ncID )
    if ( ncStatus .ne. nf_noErr ) then
        print*,nf_strError( ncStatus )
        stop
    endif

enddo


!================================================
return
stop
end subroutine outputAnalysis
