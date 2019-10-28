
include 'systemUtility/sub_build_pressureU.f90'
include 'systemUtility/sub_build_pressureV.f90'
include 'systemUtility/sub_convertGPHToGMH.f90'
include 'systemUtility/sub_build_meanDomain.f90'
include 'systemUtility/sub_check_ifObsInsideHorizontalDomain.f90'
include 'systemUtility/sub_check_ifObsInsideVerticalDomain.f90'
include 'systemUtility/sub_turnObsWithInvalidValueIntoUnavailable.f90'
include 'systemUtility/sub_rearrangeObsToBeBottomToTop.f90'
include 'systemUtility/sub_initializeAnalysis.f90'
include 'systemUtility/sub_mapObsToEachMassGrid.f90'
include 'systemUtility/sub_mapObsToEachUGrid.f90'
include 'systemUtility/sub_mapObsToEachVGrid.f90'
include 'systemUtility/sub_mapObsToEachWGrid.f90'
include 'systemUtility/sub_convertBackgroundToTemperature.f90'
include 'systemUtility/sub_convertBackgroundToRH.f90'
include 'systemUtility/sub_convertBackgroundToRefractivity.f90'
include 'systemUtility/sub_convertBackgroundToSounding.f90'
include 'systemUtility/sub_convertBackgroundToAirep.f90'
include 'systemUtility/sub_convertBackgroundToAMV.f90'
include 'systemUtility/sub_convertBackgroundToGPSRO.f90'
include 'systemUtility/sub_convertBackgroundToAIRS.f90'
include 'systemUtility/sub_convertBackgroundToQuikSCAT.f90'
include 'systemUtility/sub_convertBackgroundToASCAT.f90'
include 'systemUtility/sub_convertBackgroundToIASI.f90'
include 'systemUtility/sub_convertBackgroundToOSCAT.f90'
include 'systemUtility/sub_convertBackgroundToWindSat.f90'
include 'systemUtility/sub_setSoundingError.f90'
include 'systemUtility/sub_setAirepError.f90'
include 'systemUtility/sub_setAMVError.f90'
include 'systemUtility/sub_setGPSROError.f90'
include 'systemUtility/sub_mapGPSROFromGMHToP.f90'
include 'systemUtility/sub_setAIRSError.f90'
include 'systemUtility/sub_setQuikSCATError.f90'
include 'systemUtility/sub_setASCATError.f90'
include 'systemUtility/sub_setIASIError.f90'
include 'systemUtility/sub_setOSCATError.f90'
include 'systemUtility/sub_setWindSatError.f90'
include 'systemUtility/sub_mergeObs.f90'
include 'systemUtility/func_errorFactor.f90'
include 'systemUtility/sub_deallocate_obsListOfEachGrid.f90'

module systemUtility


interface

    subroutine build_pressureU(domain,domainSize)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)             :: domainSize
      type(domainInfo),intent(inout) :: domain(domainSize)
    end subroutine build_pressureU

    subroutine build_pressureV(domain,domainSize)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)             :: domainSize
      type(domainInfo),intent(inout) :: domain(domainSize)
    end subroutine build_pressureV

    subroutine convertGPHtoGMH(domain,domainSize)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)             :: domainSize
      type(domainInfo),intent(inout) :: domain(domainSize)
    end subroutine convertGPHtoGMH

    subroutine build_meanDomain(domain,ensembleSize,domain_mean)
      use derivedType
      implicit none
      integer,intent(in)          :: ensembleSize
      type(domainInfo),intent(in) :: domain(ensembleSize)
      type(domainInfo)            :: domain_mean
    end subroutine build_meanDomain

    subroutine check_ifObsInsideHorizontalDomain(domain,obs)
      use derivedType
      use basicUtility
      implicit none
      type(domainInfo),intent(in) :: domain
      type(obsParent),intent(inout) :: obs
    end subroutine check_ifObsInsideHorizontalDomain

    subroutine check_ifObsInsideVerticalDomain(domain,domainSize,obs)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in) :: domainSize
      type(domainInfo),intent(in) :: domain(domainSize)
      type(obsParent),intent(inout) :: obs
    end subroutine check_ifObsInsideVerticalDomain

    subroutine turnObsWithInvalidValueIntoUnavailable(obs)
      use derivedType
      implicit none
      type(obsParent),intent(inout) :: obs
    end subroutine turnObsWithInvalidValueIntoUnavailable

    subroutine rearrangeObsToBeBottomToTop(obs)
      use derivedType
      implicit none
      type(obsParent),intent(inout) :: obs
    end subroutine rearrangeObsToBeBottomToTop

    subroutine initializeAnalysis(background,analysis,ensembleSize,domain)
      use derivedType
      implicit none
      integer,intent(in) :: ensembleSize
      type(backgroundInfo),intent(in)    :: background(ensembleSize)
      type(backgroundInfo),intent(inout) :: analysis(ensembleSize)
      type(domainInfo),intent(in)        :: domain(ensembleSize)
    end subroutine initializeAnalysis

    subroutine mapObsToEachMassGrid(obsListOfEachGrid,obs,domain,systemParameters)
      use derivedType
      use basicUtility
      implicit none
      type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid  ! intent(out)
      type(obsParent),intent(inout) :: obs
      type(domainInfo),intent(in) :: domain
      type(systemParameter),intent(in) :: systemParameters
    end subroutine mapObsToEachMassGrid

    subroutine mapObsToEachUGrid(obsListOfEachGrid,obs,domain,systemParameters)
      use derivedType
      use basicUtility
      implicit none
      type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid  ! intent(out)
      type(obsParent),intent(inout) :: obs
      type(domainInfo),intent(in) :: domain
      type(systemParameter),intent(in) :: systemParameters
    end subroutine mapObsToEachUGrid

    subroutine mapObsToEachVGrid(obsListOfEachGrid,obs,domain,systemParameters)
      use derivedType
      use basicUtility
      implicit none
      type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid  ! intent(out)
      type(obsParent),intent(inout) :: obs
      type(domainInfo),intent(in) :: domain
      type(systemParameter),intent(in) :: systemParameters
    end subroutine mapObsToEachVGrid

    subroutine mapObsToEachWGrid(obsListOfEachGrid,obs,domain,systemParameters)
      use derivedType
      use basicUtility
      implicit none
      type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid  ! intent(out)
      type(obsParent),intent(inout) :: obs
      type(domainInfo),intent(in) :: domain
      type(systemParameter),intent(in) :: systemParameters
    end subroutine mapObsToEachWGrid

    subroutine convertBackgroundToTemperature(background,ensembleSize,domain)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
    end subroutine convertBackgroundToTemperature

    subroutine convertBackgroundToRH(background,ensembleSize,domain)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
    end subroutine convertBackgroundToRH

    subroutine convertBackgroundToRefractivity(background,ensembleSize,domain)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
    end subroutine convertBackgroundToRefractivity

    subroutine convertBackgroundToSounding(background,ensembleSize,domain,sounding)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(obsParent),intent(inout)   :: sounding
    end subroutine convertBackgroundToSounding

    subroutine convertBackgroundToAirep(background,ensembleSize,domain,airep)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(obsParent),intent(inout)   :: airep
    end subroutine convertBackgroundToAirep

    subroutine convertBackgroundToAMV(background,ensembleSize,domain,amv)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(obsParent),intent(inout)   :: amv
    end subroutine convertBackgroundToAMV

    subroutine convertBackgroundToGPSRO(background,ensembleSize,domain,gpsro)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(obsParent),intent(inout)   :: gpsro
    end subroutine convertBackgroundToGPSRO

    subroutine convertBackgroundToAIRS(background,ensembleSize,domain,airs)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(obsParent),intent(inout)   :: airs
    end subroutine convertBackgroundToAIRS

    subroutine convertBackgroundToQuikSCAT(background,ensembleSize,domain,domain_mean,quikscat)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(domainInfo),intent(in)     :: domain_mean
      type(obsParent),intent(inout)   :: quikscat
    end subroutine convertBackgroundToQuikSCAT

    subroutine convertBackgroundToASCAT(background,ensembleSize,domain,domain_mean,ascat)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(domainInfo),intent(in)     :: domain_mean
      type(obsParent),intent(inout)   :: ascat
    end subroutine convertBackgroundToASCAT

    subroutine convertBackgroundToIASI(background,ensembleSize,domain,iasi)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(obsParent),intent(inout)   :: iasi
    end subroutine convertBackgroundToIASI

    subroutine convertBackgroundToOSCAT(background,ensembleSize,domain,domain_mean,oscat)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(domainInfo),intent(in)     :: domain_mean
      type(obsParent),intent(inout)   :: oscat
    end subroutine convertBackgroundToOSCAT

    subroutine convertBackgroundToWindSat(background,ensembleSize,domain,domain_mean,windsat)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(domainInfo),intent(in)     :: domain_mean
      type(obsParent),intent(inout)   :: windsat
    end subroutine convertBackgroundToWindSat

    subroutine setSoundingError(sounding)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: sounding
    end subroutine setSoundingError

    subroutine setAirepError(airep)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: airep
    end subroutine setAirepError

    subroutine setAMVError(amv)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: amv
    end subroutine setAMVError

    subroutine setGPSROError(gpsro)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: gpsro
    end subroutine setGPSROError

    subroutine mapGPSROFromGMHToP(domain_mean,gpsro)
      use derivedType
      use basicUtility
      implicit none
      type(domainInfo),intent(in)     :: domain_mean
      type(obsParent),intent(inout)   :: gpsro
    end subroutine mapGPSROFromGMHToP

    subroutine setAIRSError(airs)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: airs
    end subroutine setAIRSError

    subroutine setQuikSCATError(quikscat)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: quikscat
    end subroutine setQuikSCATError

    subroutine setASCATError(ascat)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: ascat
    end subroutine setASCATError

    subroutine setIASIError(iasi)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: iasi
    end subroutine setIASIError

    subroutine setOSCATError(oscat)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: oscat
    end subroutine setOSCATError

    subroutine setWindSatError(windsat)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: windsat
    end subroutine setWindSatError

    subroutine mergeObs(allObs,obsToBeMerged)
      use derivedType
      implicit none
      type(obsParent),intent(inout) :: allObs
      type(obsParent),intent(inout) :: obsToBeMerged
    end subroutine mergeObs

    real(kind=8) function errorFactor(rc_h,rc_z,dh,dz)
      implicit none
      real(kind=8) rc_h,rc_z
      real(kind=8) dh,dz
    end function errorFactor

    subroutine deallocate_obsListOfEachGrid(obsListOfEachGrid)
      use derivedType
      implicit none
      type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid
    end subroutine deallocate_obsListOfEachGrid

end interface


end module systemUtility
