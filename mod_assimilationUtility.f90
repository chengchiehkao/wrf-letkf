

include 'assimilationUtility/sub_assimilate.f90'
include 'assimilationUtility/sub_LETKF.f90'


module assimilationUtility


interface

    subroutine assimilate(background,analysis,ensembleSize,domain,obs,obsListOfEachGrid)
      use derivedType
      use basicUtility
      use systemUtility
      use math
      implicit none
      integer,intent(in) :: ensembleSize
      type(backgroundInfo),intent(in)  :: background(ensembleSize)
      type(backgroundInfo),intent(out) :: analysis(ensembleSize)
      type(domainInfo),intent(in)      :: domain(ensembleSize)
      type(obsParent) :: obs
      type(integerVector),pointer :: obsListOfEachGrid(:,:,:)
    end subroutine assimilate

    subroutine LETKF( xb_mean , xb_pert , &
                      xa_mean , xa_pert , &
                      yo , yb , R , &
                      m , n , k )
      use math
      implicit none
      integer,parameter :: realKind=8
      integer m,n,k  ! m=model state, n=observation number, k=ensemble size
      real(realKind),dimension(m,1) :: xb_mean,xa_mean
      real(realKind),dimension(m,k) :: xb_pert,xa_pert
      real(realKind),dimension(n,1) :: yo
      real(realKind),dimension(n,k) :: yb
      real(realKind),dimension(n,n) :: R
    end subroutine LETKF

end interface

end module assimilationUtility
