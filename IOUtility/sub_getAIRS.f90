
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getAIRS(airs,systemParameters)

use derivedType

implicit none
type(obsParent),intent(out) :: airs
type(systemParameter),intent(in) :: systemParameters

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2
character(len=100) header
real(kind=8),parameter :: invalidValue = -888888.d0

integer ioStatus,fileID
integer i,io,iObsVar  ! loop counter

integer,external :: availableFileID
!================================================

airs%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: AIRS'
open(fileID,file='input/obs_input.AIRS',status='old')
do 
    header = ''
    read(fileID,'(a100)',iostat=ioStatus) header
    if ( ioStatus .lt. 0 ) then
        exit
    endif

    read(header,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum

    airs%obsNum = airs%obsNum + systemParameters%varListSize_airs*levNum +1  ! +1 for PSFC
    do i = 1,levNum+1  ! +1 for PSFC
        read(fileID,*)  ! just read through lines which are not headers.
    enddo

enddo

! Then, allocate all necessary components.
allocate( airs%obs(airs%obsNum) )

! Initilize the components just allocated.
airs%obs(:)%lon   = 0.0d0
airs%obs(:)%lat   = 0.0d0
airs%obs(:)%z     = 0.0d0
airs%obs(:)%value = 0.0d0
airs%obs(:)%error = 0.0d0
airs%obs(:)%type       = 'AIRS      '
airs%obs(:)%instrument = 'AIRS      '
airs%obs(:)%zName      = ''
airs%obs(:)%varName    = ''
airs%obs(:)%available  = .true.
do io = 1 , airs%obsNum
    allocate( airs%obs(io)%insideHorizontalDomain(systemParameters%max_domain) )
    airs%obs(io)%insideHorizontalDomain(:) = .false.
enddo

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    header = ''
    read(fileID,'(a100)',iostat=ioStatus) header
    if ( ioStatus .lt. 0 ) exit
    read(header,*) lon,lat,levNum

    serialNum = serialNum + 1
    airs%obs(serialNum)%lon = lon
    airs%obs(serialNum)%lat = lat
    airs%obs(serialNum)%zName = 'SURFACE   '
    airs%obs(serialNum)%varName = 'PSFC      '
    read(fileID,*) airs%obs(serialNum)%value
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) airs%obs(serialNum)%z,airs%obs(serialNum:serialNum+systemParameters%varListSize_airs-1)%value

        airs%obs(serialNum:serialNum+systemParameters%varListSize_airs-1)%lon = lon
        airs%obs(serialNum:serialNum+systemParameters%varListSize_airs-1)%lat = lat
        airs%obs(serialNum:serialNum+systemParameters%varListSize_airs-1)%z   = airs%obs(serialNum)%z
        airs%obs(serialNum:serialNum+systemParameters%varListSize_airs-1)%zName = 'P         '
        airs%obs(serialNum:serialNum+systemParameters%varListSize_airs-1)%varName = systemParameters%varList_airs(1:systemParameters%varListSize_airs)
        serialNum = serialNum + (systemParameters%varListSize_airs-1)
    enddo

enddo


! Convert temperature in Celsius to Kelvin
do io = 1 , airs%obsNum
    if ( trim(adjustl(airs%obs(io)%varName)).eq.'T' .and. &
         dabs(airs%obs(io)%value - invalidValue) .gt. dabs(invalidValue*epsilon(1.d0)) ) then
        airs%obs(io)%value = airs%obs(io)%value + 273.15d0
    endif
enddo

do io = 1 , airs%obsNum
    if ( trim(adjustl(airs%obs(io)%zName)) .ne. 'SURFACE' ) then
        airs%obs(io)%z = 100.d0 * airs%obs(io)%z  !  convert hPa to Pa
    else if ( trim(adjustl(airs%obs(io)%varName)) .eq. 'PSFC' ) then
        airs%obs(io)%available = .false.  ! Disable PSFC by default because it is actually an uknown field in the retrival of AIRS.
    endif
enddo

do iObsVar = 1 , systemParameters%varListSize_airs
    where ( adjustl(airs%obs(:)%varName) .eq. adjustl(systemParameters%varList_airs(iObsVar)) )
        airs%obs%available = systemParameters%use_varList_airs(iObsVar)
    end where
enddo


close(fileID)
!================================================
return
stop
end subroutine getAIRS
