
include 'basicUtility/func_inPolygon.f90'
include 'basicUtility/func_isCollinear.f90'
include 'basicUtility/sub_interp1d.f90'
include 'basicUtility/sub_interp2d.f90'
include 'basicUtility/sub_locateAsIndex1d.f90'
include 'basicUtility/sub_locateAsIndex2d.f90'
include 'basicUtility/sub_quickSortWithIndex.f90'
include 'basicUtility/func_greatCircleDistance.f90'
include 'basicUtility/func_greatCircleDistance_preCalc.f90'
include 'basicUtility/func_convertThetaAndPToTemperature.f90'
include 'basicUtility/func_convertTAndPAndQvToRH.f90'


module basicUtility

interface

    logical function inPolygon(xVertex,yVertex,vertexSize,xTarget,yTarget)
      implicit none
      integer,intent(in)                    :: vertexSize     ! size of vertices.
      real,intent(in),dimension(vertexSize) :: xVertex,yVertex  ! x-value and y-value of vertices.
      real,intent(in)                       :: xTarget,yTarget  ! x-value and y-value of target point.
    end function inPolygon

    logical function isCollinear(xEndpoint,yEndpoint,xTarget,yTarget)
      implicit none
      real,intent(in),dimension(2) :: xEndpoint,yEndpoint  ! Endpoints of the line segment.
      real,intent(in)              :: xTarget  ,yTarget    ! target point.
    end function isCollinear

    subroutine interp1d( xRef    , yRef    , refSize   , &
                         xTarget , yTarget , targetSize, &
                         method , extrap , invalidValue )
      implicit none
      integer,intent(in)                    :: refSize    ! size of referenced values.
      real,intent(in),dimension(refSize)    :: xRef,yRef  ! referenced x & y.
      integer,intent(in)                     :: targetSize  ! size of target values.
      real,intent(in) ,dimension(targetSize) :: xTarget     ! target x.
      real,intent(out),dimension(targetSize) :: yTarget     ! target y.
      integer,intent(in) :: method  ! method=1 means nearest; method=2 means linear.
      logical,intent(in) :: extrap
      real,intent(in)    :: invalidValue
    end subroutine interp1d

    subroutine interp2d( xRef           , yRef           , zRef , &
                         refRankOneSize , refRankTwoSize , &
                         xTarget , yTarget , zTarget , targetSize, &
                         method , invalidValue )
      implicit none
      integer,intent(in)                                       :: refRankOneSize,refRankTwoSize  ! size of referenced values.
      real,intent(in),dimension(refRankOneSize,refRankTwoSize) :: xRef,yRef,zRef  ! referenced x & y.
      integer,intent(in)                     :: targetSize  ! size of target values.
      real,intent(in) ,dimension(targetSize) :: xTarget,yTarget  ! target x & y.
      real,intent(out),dimension(targetSize) :: zTarget  ! target z.
      integer,intent(in) :: method  ! 1=bilinear
      real,intent(in)    :: invalidValue
    end subroutine interp2d

    subroutine locateAsIndex1d(valueRef,valueRefSize,valueTarget,indexReturned)
      implicit none
      integer,intent(in)                      :: valueRefSize  ! size of referenced values.
      real,intent(in),dimension(valueRefSize) :: valueRef      ! referenced value array.
      real,intent(in)                         :: valueTarget   ! target value.
      integer,intent(out)                     :: indexReturned
    end subroutine locateAsIndex1d

    subroutine locateAsIndex2d(xRef,yRef,refRankOneSize,refRankTwoSize, &
                               xTarget,yTarget, &
                               rankOneIndexReturned,rankTwoIndexReturned)
      implicit none
      integer,intent(in)                      :: refRankOneSize,refRankTwoSize  ! size of 1st and 2nd rank of the grids.
      real,intent(in),dimension(refRankOneSize,refRankTwoSize) :: xRef,yRef     ! x- and y-value of the grids.
      real,intent(in)                         :: xTarget,yTarget  ! x- and y-value of the target point.
      integer,intent(out)                     :: rankOneIndexReturned,rankTwoIndexReturned
    end subroutine locateAsIndex2d

    recursive subroutine quickSortWithIndex(dataToBeSorted,sizeOfData,indexOfData)
      implicit none
      integer,intent(in)    :: sizeOfData
      real,intent(inout)    :: dataToBeSorted(sizeOfData)
      integer,intent(inout) :: indexOfData(sizeOfData)
    end subroutine quickSortWithIndex

    real function greatCircleDistance(point1,point2)
      implicit none
      real,intent(in) :: point1(2),point2(2)
    end function greatCircleDistance

    real function greatCircleDistance_preCalc( point1_lon     , point2_lon     , &
                                               sin_point1_lat , sin_point2_lat , &
                                               cos_point1_lat , cos_point2_lat )
      implicit none
      real,intent(in) :: point1_lon,point2_lon
      real,intent(in) :: sin_point1_lat,sin_point2_lat
      real,intent(in) :: cos_point1_lat,cos_point2_lat
    end function greatCircleDistance_preCalc

    real function convertThetaAndPToTemperature(theta,p)
      implicit none
      real,intent(in) :: theta  ! potential temperature (unit: K)
      real,intent(in) :: p      ! pressure (unit: Pa)
    end function convertThetaAndPToTemperature

    real function convertTAndPAndQvToRH( t_k , p , qv )
      implicit none
      real,intent(in) :: t_k  ! normal temperature (unit: K)
      real,intent(in) :: p    ! pressure (unit: Pa)
      real,intent(in) :: qv   ! water vapor mixing ratio (unit: kg/kg)
    end function convertTAndPAndQvToRH


end interface


end module basicUtility
