

include 'IOUtility/func_availableFileID.f90'
include 'IOUtility/sub_getAirep.f90'
include 'IOUtility/sub_getAIRS.f90'
include 'IOUtility/sub_getAMV.f90'
include 'IOUtility/sub_getASCAT.f90'
include 'IOUtility/sub_getBackground.f90'
include 'IOUtility/sub_getCYGNSS.f90'
include 'IOUtility/sub_getDomain.f90'
include 'IOUtility/sub_getGPSRO.f90'
include 'IOUtility/sub_getIASI.f90'
include 'IOUtility/sub_getOSCAT.f90'
include 'IOUtility/sub_getQuikSCAT.f90'
include 'IOUtility/sub_getSounding.f90'
include 'IOUtility/sub_getSynop.f90'
include 'IOUtility/sub_getSystemParameter.f90'
include 'IOUtility/sub_getWindSat.f90'
include 'IOUtility/sub_outputAnalysis.f90'

module IOUtility


interface

    integer function availableFileID()
      implicit none
    end function availableFileID

    subroutine getAirep(airep,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: airep
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getAirep

    subroutine getAIRS(airs,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: airs
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getAIRS

    subroutine getAMV(amv,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: amv
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getAMV

    subroutine getASCAT(ascat,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: ascat
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getASCAT

    subroutine getBackground(background,ensembleSize,domain,domainID)
      use derivedType
      implicit none
      include 'netcdf.inc'
      integer,intent(in)                 :: ensembleSize
      type(backgroundInfo),intent(inout) :: background(ensembleSize)
      type(domainInfo),intent(in)        :: domain(ensembleSize)
      integer,intent(in)                 :: domainID
    end subroutine getBackground

    subroutine getCYGNSS(cygnss,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: cygnss
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getCYGNSS

    subroutine getDomain(domain,domainSize,domainID)
      use derivedType
      implicit none
      include 'netcdf.inc'
      integer,intent(in)           :: domainSize
      type(domainInfo),intent(out) :: domain(domainSize)
      integer,intent(in)           :: domainID
    end subroutine getDomain

    subroutine getGPSRO(GPSRO,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: GPSRO
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getGPSRO

    subroutine getIASI(iasi,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: iasi
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getIASI

    subroutine getOSCAT(oscat,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: oscat
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getOSCAT

    subroutine getQuikSCAT(quikscat,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: quikscat
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getQuikSCAT

    subroutine getSounding(sounding,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: sounding
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getSounding

    subroutine getSynop(synop,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: synop
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getSynop
 
    subroutine getSystemParameter(systemParameters)
      use derivedType
      implicit none
      type(systemParameter),intent(out) :: systemParameters
    end subroutine getSystemParameter

    subroutine getWindSat(windsat,systemParameters)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: windsat
      type(systemParameter),intent(in) :: systemParameters
    end subroutine getWindSat

    subroutine outputAnalysis(analysis,ensembleSize,domain,domainID)
      use derivedType
      implicit none
      include 'netcdf.inc'
      integer,intent(in)                 :: ensembleSize
      type(backgroundInfo),intent(inout) :: analysis(ensembleSize)
      type(domainInfo),intent(in)        :: domain(ensembleSize)
      integer,intent(in)                 :: domainID
    end subroutine outputAnalysis

end interface

end module IOUtility
