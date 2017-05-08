

include 'IOUtility/func_availableFileID.f90'
include 'IOUtility/sub_getAirep.f90'
include 'IOUtility/sub_getAIRS.f90'
include 'IOUtility/sub_getAMV.f90'
include 'IOUtility/sub_getASCAT.f90'
include 'IOUtility/sub_getBackground.f90'
include 'IOUtility/sub_getDomain.f90'
include 'IOUtility/sub_getGPSRO.f90'
include 'IOUtility/sub_getIASI.f90'
include 'IOUtility/sub_getOSCAT.f90'
include 'IOUtility/sub_getQuikSCAT.f90'
include 'IOUtility/sub_getSounding.f90'
include 'IOUtility/sub_getSynop.f90'
include 'IOUtility/sub_getSystemParameter.f90'
include 'IOUtility/sub_outputAnalysis.f90'

module IOUtility


interface

    integer function availableFileID()
      implicit none
    end function availableFileID

    subroutine getAirep(airep,varList,varListSize,use_varList)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: airep
      integer,intent(in)                                   :: varListSize
      character(len=10),dimension(varListSize),intent(in)  :: varList
      logical,dimension(varListSize),intent(in)            :: use_varList
    end subroutine getAirep

    subroutine getAIRS(airs,varList,varListSize,use_varList)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: airs
      integer,intent(in)                                   :: varListSize
      character(len=10),dimension(varListSize),intent(in)  :: varList
      logical,dimension(varListSize),intent(in)            :: use_varList
    end subroutine getAIRS

    subroutine getAMV(amv,varList,varListSize,use_varList)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: amv
      integer,intent(in)                                   :: varListSize
      character(len=10),dimension(varListSize),intent(in)  :: varList
      logical,dimension(varListSize),intent(in)            :: use_varList
    end subroutine getAMV

    subroutine getASCAT(ascat,varList,varListSize,use_varList)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: ascat
      integer,intent(in)                                   :: varListSize
      character(len=10),dimension(varListSize),intent(in)  :: varList
      logical,dimension(varListSize),intent(in)            :: use_varList
    end subroutine getASCAT

    subroutine getBackground(background,ensembleSize,domain)
      use derivedType
      implicit none
      include 'netcdf.inc'
      integer,intent(in)                 :: ensembleSize
      type(backgroundInfo),intent(inout) :: background(ensembleSize)
      type(domainInfo),intent(in)        :: domain(ensembleSize)
    end subroutine getBackground

    subroutine getDomain(domain,domainSize)
      use derivedType
      implicit none
      include 'netcdf.inc'
      integer,intent(in)           :: domainSize
      type(domainInfo),intent(out) :: domain(domainSize)
    end subroutine getDomain

    subroutine getGPSRO(GPSRO,varList,varListSize,use_varList)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: GPSRO
      integer,intent(in)                                   :: varListSize
      character(len=10),dimension(varListSize),intent(in)  :: varList
      logical,dimension(varListSize),intent(in)            :: use_varList
    end subroutine getGPSRO

    subroutine getIASI(iasi,varList,varListSize,use_varList)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: iasi
      integer,intent(in)                                   :: varListSize
      character(len=10),dimension(varListSize),intent(in)  :: varList
      logical,dimension(varListSize),intent(in)            :: use_varList
    end subroutine getIASI

    subroutine getOSCAT(oscat,varList,varListSize,use_varList)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: oscat
      integer,intent(in)                                   :: varListSize
      character(len=10),dimension(varListSize),intent(in)  :: varList
      logical,dimension(varListSize),intent(in)            :: use_varList
    end subroutine getOSCAT

    subroutine getQuikSCAT(quikscat,varList,varListSize,use_varList)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: quikscat
      integer,intent(in)                                   :: varListSize
      character(len=10),dimension(varListSize),intent(in)  :: varList
      logical,dimension(varListSize),intent(in)            :: use_varList
    end subroutine getQuikSCAT

    subroutine getSounding(sounding,varList,varListSize,use_varList)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: sounding
      integer,intent(in)                                   :: varListSize
      character(len=10),dimension(varListSize),intent(in)  :: varList
      logical,dimension(varListSize),intent(in)            :: use_varList
    end subroutine getSounding

    subroutine getSynop(synop,varList,varListSize,use_varList)
      use derivedType
      implicit none
      type(obsParent),intent(out) :: synop
      integer,intent(in)                                   :: varListSize
      character(len=10),dimension(varListSize),intent(in)  :: varList
      logical,dimension(varListSize),intent(in)            :: use_varList
    end subroutine getSynop
 
    subroutine getSystemParameter(systemParameters)
      use derivedType
      implicit none
      type(systemParameter),intent(out) :: systemParameters
    end subroutine getSystemParameter

    subroutine outputAnalysis(analysis,ensembleSize,domain)
      use derivedType
      implicit none
      include 'netcdf.inc'
      integer,intent(in)                 :: ensembleSize
      type(backgroundInfo),intent(inout) :: analysis(ensembleSize)
      type(domainInfo),intent(in)        :: domain(ensembleSize)
    end subroutine outputAnalysis

end interface

end module IOUtility
