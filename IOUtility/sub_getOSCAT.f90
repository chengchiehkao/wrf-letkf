
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getOSCAT(oscat,systemParameters)

use derivedType

implicit none
type(obsParent),intent(out) :: oscat
type(systemParameter),intent(in) :: systemParameters

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2

integer ioStatus,fileID
integer i,io,iObsVar  ! loop counter

integer,external :: availableFileID
!================================================

oscat%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: OSCAT'
open(fileID,file='input/obs_input.oscat',status='old')
do 
    read(fileID,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum
    if ( ioStatus .lt. 0 ) then
        exit
    endif
    oscat%obsNum = oscat%obsNum + systemParameters%varListSize_oscat*levNum
    do i = 1,levNum
        read(fileID,*)  ! just read through lines which are not headers.
    enddo
enddo

! Then, allocate all necessary components.
allocate( oscat%obs(oscat%obsNum) )

! Initilize the components just allocated.
oscat%obs(:)%lon   = 0.0d0
oscat%obs(:)%lat   = 0.0d0
oscat%obs(:)%z     = 0.0d0
oscat%obs(:)%value = 0.0d0
oscat%obs(:)%error = 0.0d0
oscat%obs(:)%instrument = 'OSCAT'
oscat%obs(:)%zName      = 'P'
oscat%obs(:)%varName    = ''
oscat%obs(:)%available  = .true.
do io = 1 , oscat%obsNum
    allocate( oscat%obs(io)%insideHorizontalDomain(systemParameters%max_domain) )
    oscat%obs(io)%insideHorizontalDomain(:) = .false.
enddo

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) oscat%obs(serialNum)%z,oscat%obs(serialNum:serialNum+systemParameters%varListSize_oscat-1)%value

        oscat%obs(serialNum:serialNum+systemParameters%varListSize_oscat-1)%lon = lon
        oscat%obs(serialNum:serialNum+systemParameters%varListSize_oscat-1)%lat = lat
        oscat%obs(serialNum:serialNum+systemParameters%varListSize_oscat-1)%z   = oscat%obs(serialNum)%z
        oscat%obs(serialNum:serialNum+systemParameters%varListSize_oscat-1)%varName = systemParameters%varList_oscat(1:systemParameters%varListSize_oscat)
        serialNum = serialNum + (systemParameters%varListSize_oscat-1)
    enddo
enddo


oscat%obs(:)%z = 100.d0 * oscat%obs(:)%z  !  convert hPa to Pa

do iObsVar = 1 , systemParameters%varListSize_oscat
    where ( adjustl(oscat%obs(:)%varName) .eq. adjustl(systemParameters%varList_oscat(iObsVar)) )
        oscat%obs%available = systemParameters%use_varList_oscat(iObsVar)
    end where
enddo


close(fileID)
!================================================
return
stop
end subroutine getOSCAT
