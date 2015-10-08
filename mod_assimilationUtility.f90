

include 'assimilationUtility/sub_assimilate_massGrid.f90'
include 'assimilationUtility/sub_assimilate_uGrid.f90'
include 'assimilationUtility/sub_assimilate_vGrid.f90'
include 'assimilationUtility/sub_assimilate_wGrid.f90'
include 'assimilationUtility/sub_LETKF.f90'


module assimilationUtility


interface

    subroutine assimilate_massGrid(background,analysis,ensembleSize,domain,domain_mean,obs,obsListOfEachGrid,systemParameters)
      use derivedType
      use basicUtility
      use systemUtility
      use math
      implicit none
      integer,intent(in) :: ensembleSize
      type(backgroundInfo),intent(in)  :: background(ensembleSize)
      type(backgroundInfo),intent(inout) :: analysis(ensembleSize)
      type(domainInfo),intent(in)      :: domain(ensembleSize)
      type(domainInfo),intent(in)      :: domain_mean
      type(obsParent) :: obs
      type(integerVector),pointer :: obsListOfEachGrid(:,:,:)
      type(systemParameter),intent(in) :: systemParameters
    end subroutine assimilate_massGrid

    subroutine assimilate_uGrid(background,analysis,ensembleSize,domain,domain_mean,obs,obsListOfEachGrid,systemParameters)
      use derivedType
      use basicUtility
      use systemUtility
      use math
      implicit none
      integer,intent(in) :: ensembleSize
      type(backgroundInfo),intent(in)  :: background(ensembleSize)
      type(backgroundInfo),intent(inout) :: analysis(ensembleSize)
      type(domainInfo),intent(in)      :: domain(ensembleSize)
      type(domainInfo),intent(in)      :: domain_mean
      type(obsParent) :: obs
      type(integerVector),pointer :: obsListOfEachGrid(:,:,:)
      type(systemParameter),intent(in) :: systemParameters
    end subroutine assimilate_uGrid

    subroutine assimilate_vGrid(background,analysis,ensembleSize,domain,domain_mean,obs,obsListOfEachGrid,systemParameters)
      use derivedType
      use basicUtility
      use systemUtility
      use math
      implicit none
      integer,intent(in) :: ensembleSize
      type(backgroundInfo),intent(in)  :: background(ensembleSize)
      type(backgroundInfo),intent(inout) :: analysis(ensembleSize)
      type(domainInfo),intent(in)      :: domain(ensembleSize)
      type(domainInfo),intent(in)      :: domain_mean
      type(obsParent) :: obs
      type(integerVector),pointer :: obsListOfEachGrid(:,:,:)
      type(systemParameter),intent(in) :: systemParameters
    end subroutine assimilate_vGrid

    subroutine assimilate_wGrid(background,analysis,ensembleSize,domain,domain_mean,obs,obsListOfEachGrid,systemParameters)
      use derivedType
      use basicUtility
      use systemUtility
      use math
      implicit none
      integer,intent(in) :: ensembleSize
      type(backgroundInfo),intent(in)  :: background(ensembleSize)
      type(backgroundInfo),intent(inout) :: analysis(ensembleSize)
      type(domainInfo),intent(in)      :: domain(ensembleSize)
      type(domainInfo),intent(in)      :: domain_mean
      type(obsParent) :: obs
      type(integerVector),pointer :: obsListOfEachGrid(:,:,:)
      type(systemParameter),intent(in) :: systemParameters
    end subroutine assimilate_wGrid

    subroutine LETKF( xb_mean , xb_pert , &
                      xa_mean , xa_pert , &
                      yo , yb , R , &
                      m , n , k , inflationFactor )
      use math
      implicit none
      integer,parameter :: realKind=8
      integer m,n,k  ! m=model state, n=observation number, k=ensemble size
      real(realKind),dimension(m,1) :: xb_mean,xa_mean
      real(realKind),dimension(m,k) :: xb_pert,xa_pert
      real(realKind),dimension(n,1) :: yo
      real(realKind),dimension(n,k) :: yb
      real(realKind),dimension(n)   :: R
      real(realKind) :: inflationFactor
    end subroutine LETKF

end interface

end module assimilationUtility
