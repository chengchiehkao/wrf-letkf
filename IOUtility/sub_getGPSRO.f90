
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getGPSRO(GPSRO,varList,varListSize)

use derivedType

implicit none
type(obsParent),intent(out) :: GPSRO
integer,intent(in)                                   :: varListSize
character(len=10),dimension(varListSize),intent(in)  :: varList

real(kind=8) lon,lat,rfict
integer levNum,serialNum
real :: dummyArg(1:3)

integer ioStatus,fileID
integer i  ! loop counter

integer,external :: availableFileID
!================================================

GPSRO%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: GPSRO'
open(fileID,file='input/obs_input.cosmic',status='old')
do 
    read(fileID,*,iostat=ioStatus) dummyArg(1:3) , levNum
    if ( ioStatus .lt. 0 ) then
        exit
    endif
    GPSRO%obsNum = GPSRO%obsNum + varListSize*levNum
    do i = 1,levNum
        read(fileID,*)  ! just read through lines which are not headers.
    enddo
enddo

! Then, allocate all necessary components.
allocate( GPSRO%obs(GPSRO%obsNum) )

! Initilize the components just allocated.
GPSRO%obs(:)%lon   = 0.0d0
GPSRO%obs(:)%lat   = 0.0d0
GPSRO%obs(:)%z     = 0.0d0
GPSRO%obs(:)%value = 0.0d0
GPSRO%obs(:)%error = 0.0d0
GPSRO%obs(:)%rfict = 0.0d0
GPSRO%obs(:)%instrument = 'COSMIC'
GPSRO%obs(:)%zName      = 'GMH'
GPSRO%obs(:)%varName    = ''
GPSRO%obs(:)%available  = .true.

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,rfict,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) GPSRO%obs(serialNum)%z,GPSRO%obs(serialNum:serialNum+varListSize-1)%value

        GPSRO%obs(serialNum:serialNum+varListSize-1)%lon = lon
        GPSRO%obs(serialNum:serialNum+varListSize-1)%lat = lat
        GPSRO%obs(serialNum:serialNum+varListSize-1)%rfict = rfict
        GPSRO%obs(serialNum:serialNum+varListSize-1)%z = GPSRO%obs(serialNum)%z
        GPSRO%obs(serialNum:serialNum+varListSize-1)%varName = varList(1:varListSize)
        serialNum = serialNum + (varListSize-1)
    enddo
enddo

close(fileID)
!================================================
return
stop
end subroutine getGPSRO
