
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getWindSat(windsat,systemParameters)

use derivedType

implicit none
type(obsParent),intent(out) :: windsat
type(systemParameter),intent(in) :: systemParameters

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2

integer ioStatus,fileID
integer i,io,iObsVar  ! loop counter

integer,external :: availableFileID
!================================================

windsat%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: WindSat'
open(fileID,file='input/obs_input.windsat',status='old')
do 
    read(fileID,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum
    if ( ioStatus .lt. 0 ) then
        exit
    endif
    windsat%obsNum = windsat%obsNum + systemParameters%varListSize_windsat*levNum
    do i = 1,levNum
        read(fileID,*)  ! just read through lines which are not headers.
    enddo
enddo

! Then, allocate all necessary components.
allocate( windsat%obs(windsat%obsNum) )

! Initilize the components just allocated.
windsat%obs(:)%lon   = 0.0d0
windsat%obs(:)%lat   = 0.0d0
windsat%obs(:)%z     = 0.0d0
windsat%obs(:)%value = 0.0d0
windsat%obs(:)%error = 0.0d0
windsat%obs(:)%instrument = 'WindSat'
windsat%obs(:)%zName      = 'P'
windsat%obs(:)%varName    = ''
windsat%obs(:)%available  = .true.
do io = 1 , windsat%obsNum
    allocate( windsat%obs(io)%insideHorizontalDomain(systemParameters%max_domain) )
    windsat%obs(io)%insideHorizontalDomain(:) = .false.
enddo

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) windsat%obs(serialNum)%z,windsat%obs(serialNum:serialNum+systemParameters%varListSize_windsat-1)%value

        windsat%obs(serialNum:serialNum+systemParameters%varListSize_windsat-1)%lon = lon
        windsat%obs(serialNum:serialNum+systemParameters%varListSize_windsat-1)%lat = lat
        windsat%obs(serialNum:serialNum+systemParameters%varListSize_windsat-1)%z   = windsat%obs(serialNum)%z
        windsat%obs(serialNum:serialNum+systemParameters%varListSize_windsat-1)%varName = systemParameters%varList_windsat(1:systemParameters%varListSize_windsat)
        serialNum = serialNum + (systemParameters%varListSize_windsat-1)
    enddo
enddo


windsat%obs(:)%z = 100.d0 * windsat%obs(:)%z  !  convert hPa to Pa

do iObsVar = 1 , systemParameters%varListSize_windsat
    where ( adjustl(windsat%obs(:)%varName) .eq. adjustl(systemParameters%varList_windsat(iObsVar)) )
        windsat%obs%available = systemParameters%use_varList_windsat(iObsVar)
    end where
enddo


close(fileID)
!================================================
return
stop
end subroutine getWindSat
