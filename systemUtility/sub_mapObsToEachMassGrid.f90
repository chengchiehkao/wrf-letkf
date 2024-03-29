
subroutine mapObsToEachMassGrid(obsListOfEachGrid,obs,domain,systemParameters)

use derivedType
use basicUtility

implicit none
type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid  ! intent(out)
type(obsParent),intent(inout) :: obs
type(domainInfo),intent(in) :: domain
type(systemParameter),intent(in) :: systemParameters

integer,allocatable,dimension(:) :: indexOfObs
integer :: indexBuffer(obs%obsNum)

real(kind=8) :: rd,rd_z
real(kind=8) zDiff

integer :: blockNum_rank1=36,blockNum_rank2=30,subDomainObsNum
integer,allocatable,dimension(:) :: subDomainStart_rank1,subDomainStart_rank2
integer,allocatable,dimension(:) :: subDomainEnd_rank1  ,subDomainEnd_rank2
integer,allocatable,dimension(:) :: subDomainIndexBuffer
real(kind=8) :: subDomainCenter_lon,subDomainCenter_lat , longestDistanceBetweenCenterAndGrids
real(kind=8) :: subDomainCenter_lonInRad , sin_subDomainCenter_latInRad , cos_subDomainCenter_latInRad

real(kind=8),allocatable :: zGrid_preprocessed(:,:,:),zObs_preprocessed(:)

real(kind=8),allocatable :: domainLonInRad(:,:) , obsLonInRad(:)
real(kind=8),allocatable :: sin_domainLatInRad(:,:) , sin_obsLatInRad(:)
real(kind=8),allocatable :: cos_domainLatInRad(:,:) , cos_obsLatInRad(:)
real(kind=8) :: zLayerMax,zLayerMin

real(kind=8),parameter :: radian = 57.29577951308232d0

integer io,iwe,isn,iz,ir1,ir2
!================================================

rd    = systemParameters % rd
rd_z  = systemParameters % rd_z


allocate( subDomainIndexBuffer( obs%obsNum ) )
allocate( subDomainStart_rank1( blockNum_rank1 ) , subDomainStart_rank2( blockNum_rank2 ) )
allocate( subDomainEnd_rank1( blockNum_rank1 )   , subDomainEnd_rank2( blockNum_rank2 ) )
allocate( obsListOfEachGrid( domain%size_westToEast , domain%size_southToNorth , domain%size_bottomToTop ) )

subDomainIndexBuffer(:) = 0
subDomainStart_rank1(:) = 0
subDomainStart_rank2(:) = 0
subDomainEnd_rank1(:) = 0
subDomainEnd_rank2(:) = 0
obsListOfEachGrid(:,:,:)%vectorSize = 0


do ir1 = 1 , blockNum_rank1
    subDomainStart_rank1( ir1 ) = 1 + nint( real((ir1-1)*(domain%size_westToEast-1)) / real(blockNum_rank1) )
enddo
subDomainEnd_rank1( 1:blockNum_rank1-1 ) = subDomainStart_rank1(2:blockNum_rank1) - 1
subDomainEnd_rank1( blockNum_rank1 )     = domain%size_westToEast

do ir2 = 1 , blockNum_rank2
    subDomainStart_rank2( ir2 ) = 1 + nint( real((ir2-1)*(domain%size_southToNorth-1)) / real(blockNum_rank2) )
enddo
subDomainEnd_rank2( 1:blockNum_rank2-1 ) = subDomainStart_rank2(2:blockNum_rank2) - 1
subDomainEnd_rank2( blockNum_rank2 )     = domain%size_southToNorth



allocate(zGrid_preprocessed(domain%size_westToEast,domain%size_southToNorth,domain%size_bottomToTop))
allocate(zObs_preprocessed(obs%obsNum))

zGrid_preprocessed(:,:,:) = dlog(domain%pressure(:,:,:) * 0.01d0)
zObs_preprocessed(:)      = dlog(obs%obs(:)%z * 0.01d0)

allocate( domainLonInRad( domain%size_westToEast , domain%size_southToNorth ) , obsLonInRad( obs%obsNum ) )
allocate( sin_domainLatInRad( domain%size_westToEast , domain%size_southToNorth ) , sin_obsLatInRad( obs%obsNum ) )
allocate( cos_domainLatInRad( domain%size_westToEast , domain%size_southToNorth ) , cos_obsLatInRad( obs%obsNum ) )

domainLonInRad(:,:) = domain%lon(:,:) / radian
obsLonInRad(:)    = obs%obs(:)%lon / radian

sin_domainLatInRad(:,:) = dsin( domain%lat(:,:) / radian )
cos_domainLatInRad(:,:) = dcos( domain%lat(:,:) / radian )
sin_obsLatInRad(:)    = dsin( obs%obs(:)%lat / radian )
cos_obsLatInRad(:)    = dcos( obs%obs(:)%lat / radian )


do ir2 = 1 , blockNum_rank2
do ir1 = 1 , blockNum_rank1

    ! find out the center of sub-domain
    subDomainCenter_lon = sum( domain%lon( subDomainStart_rank1(ir1):subDomainEnd_rank1(ir1) , &
                                           subDomainStart_rank2(ir2):subDomainEnd_rank2(ir2) ) ) &
                               / real( (subDomainEnd_rank1(ir1)-subDomainStart_rank1(ir1)+1) * &
                                       (subDomainEnd_rank2(ir2)-subDomainStart_rank2(ir2)+1) )

    subDomainCenter_lat = sum( domain%lat( subDomainStart_rank1(ir1):subDomainEnd_rank1(ir1) , &
                                           subDomainStart_rank2(ir2):subDomainEnd_rank2(ir2) ) ) &
                               / real( (subDomainEnd_rank1(ir1)-subDomainStart_rank1(ir1)+1) * &
                                       (subDomainEnd_rank2(ir2)-subDomainStart_rank2(ir2)+1) )

    subDomainCenter_lonInRad = subDomainCenter_lon / radian
    sin_subDomainCenter_latInRad = dsin( subDomainCenter_lat / radian )
    cos_subDomainCenter_latInRad = dcos( subDomainCenter_lat / radian )
    

    ! determine the longest distance between sub-domain center and sub-domain grid points
    longestDistanceBetweenCenterAndGrids = -9.d15
    do isn = subDomainStart_rank2(ir2) , subDomainEnd_rank2(ir2)
    do iwe = subDomainStart_rank1(ir1) , subDomainEnd_rank1(ir1)
        longestDistanceBetweenCenterAndGrids = max( greatCircleDistance( (/ subDomainCenter_lon , subDomainCenter_lat /) , &
                                                                         (/ domain%lon(iwe,isn) , domain%lat(iwe,isn) /) ) ,&
                                                    longestDistanceBetweenCenterAndGrids )
    enddo
    enddo

    ! map obs' in sub-domain index buffer to each grids in sub-domain
!$omp parallel do default(shared) private(iwe,isn,iz,zDiff,io,indexBuffer,subDomainObsNum,subDomainIndexBuffer,zLayerMax,zLayerMin) schedule(dynamic,1)
    do iz  = 1 , domain % size_bottomToTop

        ! gather obs' into sub-domain index buffer
        subDomainObsNum = 0
        subDomainIndexBuffer(:) = 0
        zLayerMax = maxval(zGrid_preprocessed(subDomainStart_rank1(ir1):subDomainEnd_rank1(ir1),subDomainStart_rank2(ir2):subDomainEnd_rank2(ir2),iz))
        zLayerMin = minval(zGrid_preprocessed(subDomainStart_rank1(ir1):subDomainEnd_rank1(ir1),subDomainStart_rank2(ir2):subDomainEnd_rank2(ir2),iz))
        do io  = 1 , obs % obsNum
            if ( obs%obs(io)%available ) then
                if ( zObs_preprocessed(io)-zLayerMax .le. rd_z .and. zObs_preprocessed(io)-zLayerMin .ge. -(rd_z) ) then
                    if ( greatCircleDistance_preCalc( subDomainCenter_lonInRad , obsLonInRad(io) , &
                                                      sin_subDomainCenter_latInRad , sin_obsLatInRad(io) , &
                                                      cos_subDomainCenter_latInRad , cos_obsLatInRad(io) ) .le. (rd+longestDistanceBetweenCenterAndGrids) ) then
                        subDomainObsNum = subDomainObsNum + 1
                        subDomainIndexBuffer(subDomainObsNum) = io
                    endif
                endif
            endif
        enddo

        indexBuffer(:) = 0
        do isn = subDomainStart_rank2(ir2) , subDomainEnd_rank2(ir2)
        do iwe = subDomainStart_rank1(ir1) , subDomainEnd_rank1(ir1)

            zDiff = 0.d0  ! shall not less than -rd_z
            do io  = 1 , subDomainObsNum
                if ( zDiff .lt. -rd_z ) exit
                if ( obs%obs(subDomainIndexBuffer(io))%available ) then
                    zDiff = zObs_preprocessed(subDomainIndexBuffer(io)) - zGrid_preprocessed(iwe,isn,iz)
                    if ( dabs( zDiff ) .le. rd_z ) then
                        if ( greatCircleDistance_preCalc( domainLonInRad(iwe,isn) , obsLonInRad(subDomainIndexBuffer(io)) , &
                                                          sin_domainLatInRad(iwe,isn) , sin_obsLatInRad(subDomainIndexBuffer(io)) , &
                                                          cos_domainLatInRad(iwe,isn) , cos_obsLatInRad(subDomainIndexBuffer(io)) ) .le. rd ) then
                            obsListOfEachGrid(iwe,isn,iz)%vectorSize = obsListOfEachGrid(iwe,isn,iz)%vectorSize + 1
                            indexBuffer(obsListOfEachGrid(iwe,isn,iz)%vectorSize) = subDomainIndexBuffer(io)
                        endif
                    endif
                endif
            enddo

            if ( obsListOfEachGrid(iwe,isn,iz)%vectorSize .gt. 0 ) then
                allocate( obsListOfEachGrid(iwe,isn,iz)%vector(obsListOfEachGrid(iwe,isn,iz)%vectorSize) )
                obsListOfEachGrid(iwe,isn,iz)%vector(:) = indexBuffer(1:obsListOfEachGrid(iwe,isn,iz)%vectorSize)
            endif

        enddo
        enddo

    enddo
!$omp end parallel do


enddo
enddo


!================================================
return
stop
end subroutine mapObsToEachMassGrid
