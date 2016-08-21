
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getAirep(airep,varList,varListSize,use_varList)

use derivedType

implicit none
type(obsParent),intent(out) :: airep
integer,intent(in)                                   :: varListSize
character(len=10),dimension(varListSize),intent(in)  :: varList
logical,dimension(varListSize),intent(in)            :: use_varList

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2
character(len=100) header
real(kind=8),parameter :: invalidValue = -888888.d0

integer ioStatus,fileID
integer i,io,iObsVar  ! loop counter

integer,external :: availableFileID
!================================================

airep%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: airep'
open(fileID,file='input/obs_input.airep',status='old')
do 
    header = ''
    read(fileID,'(a100)',iostat=ioStatus) header
    if ( ioStatus .lt. 0 ) then
        exit
    endif

    read(header,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum
    airep%obsNum = airep%obsNum + varListSize*levNum
    do i = 1,levNum
        read(fileID,*)  ! just read through lines which are not headers.
    enddo

enddo

! Then, allocate all necessary components.
allocate( airep%obs(airep%obsNum) )

! Initilize the components just allocated.
airep%obs(:)%lon   = 0.0d0
airep%obs(:)%lat   = 0.0d0
airep%obs(:)%z     = 0.0d0
airep%obs(:)%value = 0.0d0
airep%obs(:)%error = 0.0d0
airep%obs(:)%type       = ''
airep%obs(:)%instrument = 'AIREP     '
airep%obs(:)%zName      = ''
airep%obs(:)%varName    = ''
airep%obs(:)%available  = .true.

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    header = ''
    read(fileID,'(a100)',iostat=ioStatus) header
    if ( ioStatus .lt. 0 ) exit
    read(header,*) lon,lat,levNum

    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) airep%obs(serialNum)%z,airep%obs(serialNum:serialNum+varListSize-1)%value

        airep%obs(serialNum:serialNum+varListSize-1)%instrument = 'AIREP     '
        airep%obs(serialNum:serialNum+varListSize-1)%lon = lon
        airep%obs(serialNum:serialNum+varListSize-1)%lat = lat
        airep%obs(serialNum:serialNum+varListSize-1)%z   = airep%obs(serialNum)%z
        airep%obs(serialNum:serialNum+varListSize-1)%zName = 'P         '
        airep%obs(serialNum:serialNum+varListSize-1)%varName = varList(1:varListSize)
        serialNum = serialNum + (varListSize-1)
    enddo

enddo


! Convert temperature in Celsius to Kelvin
do io = 1 , airep%obsNum
    if ( trim(adjustl(airep%obs(io)%varName)).eq.'T' .and. &
         dabs(airep%obs(io)%value - invalidValue) .gt. dabs(invalidValue*epsilon(1.d0)) ) then
        airep%obs(io)%value = airep%obs(io)%value + 273.15d0
    endif
enddo

do io = 1 , airep%obsNum
    airep%obs(io)%z = 100.d0 * airep%obs(io)%z  !  convert hPa to Pa
enddo

do iObsVar = 1 , varListSize
    where ( adjustl(airep%obs(:)%varName) .eq. adjustl(varList(iObsVar)) )
        airep%obs%available = use_varList(iObsVar)
    end where
enddo


close(fileID)
!================================================
return
stop
end subroutine getAirep
