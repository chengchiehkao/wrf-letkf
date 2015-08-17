
include 'systemUtility/sub_build_pressureU.f90'
include 'systemUtility/sub_build_pressureV.f90'
include 'systemUtility/sub_build_meanDomain.f90'
include 'systemUtility/sub_check_ifObsInsideHorizontalDomain.f90'
include 'systemUtility/sub_check_ifObsInsideVerticalDomain.f90'
include 'systemUtility/sub_turnObsWithInvalidValueIntoUnavailable.f90'
include 'systemUtility/sub_mapObsToEachMassGrid.f90'
include 'systemUtility/sub_mapObsToEachUGrid.f90'
include 'systemUtility/sub_mapObsToEachVGrid.f90'
include 'systemUtility/sub_mapObsToEachWGrid.f90'
include 'systemUtility/sub_convertBackgroundToSounding.f90'
include 'systemUtility/sub_convertBackgroundToTemperature.f90'
include 'systemUtility/sub_setSoundingError.f90'
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

    subroutine mapObsToEachMassGrid(obsListOfEachGrid,obs,domain)
      use derivedType
      use basicUtility
      implicit none
      type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid  ! intent(out)
      type(obsParent),intent(inout) :: obs
      type(domainInfo),intent(in) :: domain
    end subroutine mapObsToEachMassGrid

    subroutine mapObsToEachUGrid(obsListOfEachGrid,obs,domain)
      use derivedType
      use basicUtility
      implicit none
      type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid  ! intent(out)
      type(obsParent),intent(inout) :: obs
      type(domainInfo),intent(in) :: domain
    end subroutine mapObsToEachUGrid

    subroutine mapObsToEachVGrid(obsListOfEachGrid,obs,domain)
      use derivedType
      use basicUtility
      implicit none
      type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid  ! intent(out)
      type(obsParent),intent(inout) :: obs
      type(domainInfo),intent(in) :: domain
    end subroutine mapObsToEachVGrid

    subroutine mapObsToEachWGrid(obsListOfEachGrid,obs,domain)
      use derivedType
      use basicUtility
      implicit none
      type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid  ! intent(out)
      type(obsParent),intent(inout) :: obs
      type(domainInfo),intent(in) :: domain
    end subroutine mapObsToEachWGrid

    subroutine convertBackgroundToTemperature(background,ensembleSize,domain)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
    end subroutine convertBackgroundToTemperature

    subroutine convertBackgroundToSounding(background,ensembleSize,domain,sounding)
      use derivedType
      use basicUtility
      implicit none
      integer,intent(in)              :: ensembleSize
      type(backgroundInfo),intent(in) :: background(ensembleSize)
      type(domainInfo),intent(in)     :: domain(ensembleSize)
      type(obsParent),intent(inout)   :: sounding
    end subroutine convertBackgroundToSounding

    subroutine setSoundingError(sounding)
      use derivedType
      use basicUtility
      implicit none
      type(obsParent),intent(inout)   :: sounding
    end subroutine setSoundingError

    real(kind=8) function errorFactor(rc,rlev,dh,dz)
      implicit none
      real(kind=8) rc,rlev
      real(kind=8) dh,dz
    end function errorFactor

    subroutine deallocate_obsListOfEachGrid(obsListOfEachGrid)
      use derivedType
      implicit none
      type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid
    end subroutine deallocate_obsListOfEachGrid

end interface


end module systemUtility
