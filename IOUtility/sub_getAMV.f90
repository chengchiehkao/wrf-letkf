
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getAMV(amv,varList,varListSize)

use derivedType

implicit none
type(obsParent),intent(out) :: amv
integer,intent(in)                                   :: varListSize
character(len=10),dimension(varListSize),intent(in)  :: varList

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2

integer ioStatus,fileID
integer i  ! loop counter

integer,external :: availableFileID
!================================================

amv%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: AMV'
open(fileID,file='input/obs_input.AMV',status='old')
do 
    read(fileID,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum
    if ( ioStatus .lt. 0 ) then
        exit
    endif
    amv%obsNum = amv%obsNum + varListSize*levNum
    do i = 1,levNum
        read(fileID,*)  ! just read through lines which are not headers.
    enddo
enddo

! Then, allocate all necessary components.
allocate( amv%obs(amv%obsNum) )

! Initilize the components just allocated.
amv%obs(:)%lon   = 0.0d0
amv%obs(:)%lat   = 0.0d0
amv%obs(:)%z     = 0.0d0
amv%obs(:)%value = 0.0d0
amv%obs(:)%error = 0.0d0
amv%obs(:)%instrument = 'AMV'
amv%obs(:)%zName      = 'P'
amv%obs(:)%varName    = ''
amv%obs(:)%available  = .true.

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) amv%obs(serialNum)%z,amv%obs(serialNum:serialNum+varListSize-1)%value

        amv%obs(serialNum:serialNum+varListSize-1)%lon = lon
        amv%obs(serialNum:serialNum+varListSize-1)%lat = lat
        amv%obs(serialNum:serialNum+varListSize-1)%z   = amv%obs(serialNum)%z
        amv%obs(serialNum:serialNum+varListSize-1)%varName = varList(1:varListSize)
        serialNum = serialNum + (varListSize-1)
    enddo
enddo


amv%obs(:)%z = 100.d0 * amv%obs(:)%z  !  convert hPa to Pa

close(fileID)
!================================================
return
stop
end subroutine getAMV
