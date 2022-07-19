
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getCYGNSS(cygnss,systemParameters)

use derivedType

implicit none
type(obsParent),intent(out) :: cygnss
type(systemParameter),intent(in) :: systemParameters

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2

integer ioStatus,fileID
integer i,io,iObsVar  ! loop counter

integer,external :: availableFileID
!================================================

cygnss%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: CYGNSS'
open(fileID,file='input/obs_input.cygnss',status='old')
do 
    read(fileID,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum
    if ( ioStatus .lt. 0 ) then
        exit
    endif
    cygnss%obsNum = cygnss%obsNum + systemParameters%varListSize_cygnss*levNum
    do i = 1,levNum
        read(fileID,*)  ! just read through lines which are not headers.
    enddo
enddo

! Then, allocate all necessary components.
allocate( cygnss%obs(cygnss%obsNum) )

! Initilize the components just allocated.
cygnss%obs(:)%lon   = 0.0d0
cygnss%obs(:)%lat   = 0.0d0
cygnss%obs(:)%z     = 0.0d0
cygnss%obs(:)%value = 0.0d0
cygnss%obs(:)%error = 0.0d0
cygnss%obs(:)%instrument = 'CYGNSS'
cygnss%obs(:)%zName      = 'P'
cygnss%obs(:)%varName    = ''
cygnss%obs(:)%available  = .true.
do io = 1 , cygnss%obsNum
    allocate( cygnss%obs(io)%insideHorizontalDomain(systemParameters%max_domain) )
    cygnss%obs(io)%insideHorizontalDomain(:) = .false.
enddo

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) cygnss%obs(serialNum)%z,cygnss%obs(serialNum:serialNum+systemParameters%varListSize_cygnss-1)%value

        cygnss%obs(serialNum:serialNum+systemParameters%varListSize_cygnss-1)%lon = lon
        cygnss%obs(serialNum:serialNum+systemParameters%varListSize_cygnss-1)%lat = lat
        cygnss%obs(serialNum:serialNum+systemParameters%varListSize_cygnss-1)%z   = cygnss%obs(serialNum)%z
        cygnss%obs(serialNum:serialNum+systemParameters%varListSize_cygnss-1)%varName = systemParameters%varList_cygnss(1:systemParameters%varListSize_cygnss)
        serialNum = serialNum + (systemParameters%varListSize_cygnss-1)
    enddo
enddo


cygnss%obs(:)%z = 100.d0 * cygnss%obs(:)%z  !  convert hPa to Pa

do iObsVar = 1 , systemParameters%varListSize_cygnss
    where ( adjustl(cygnss%obs(:)%varName) .eq. adjustl(systemParameters%varList_cygnss(iObsVar)) )
        cygnss%obs%available = systemParameters%use_varList_cygnss(iObsVar)
    end where
enddo


close(fileID)
!================================================
return
stop
end subroutine getCYGNSS
