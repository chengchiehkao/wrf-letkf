
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getGPSRO(GPSRO,systemParameters)

use derivedType

implicit none
type(obsParent),intent(out) :: GPSRO
type(systemParameter),intent(in) :: systemParameters

real(kind=8) lon,lat,rfict
integer levNum,serialNum
real :: dummyArg(1:3)

integer ioStatus,fileID
integer i,io,iObsVar  ! loop counter

integer,external :: availableFileID
!================================================

GPSRO%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: GPSRO'
open(fileID,file='input/obs_input.gpsro',status='old')
do 
    read(fileID,*,iostat=ioStatus) dummyArg(1:3) , levNum
    if ( ioStatus .lt. 0 ) then
        exit
    endif
    GPSRO%obsNum = GPSRO%obsNum + systemParameters%varListSize_gpsro*levNum
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
do io = 1 , GPSRO%obsNum
    allocate( GPSRO%obs(io)%insideHorizontalDomain(systemParameters%max_domain) )
    GPSRO%obs(io)%insideHorizontalDomain(:) = .false.
enddo

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,rfict,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) GPSRO%obs(serialNum)%z,GPSRO%obs(serialNum:serialNum+systemParameters%varListSize_gpsro-1)%value

        GPSRO%obs(serialNum:serialNum+systemParameters%varListSize_gpsro-1)%lon = lon
        GPSRO%obs(serialNum:serialNum+systemParameters%varListSize_gpsro-1)%lat = lat
        GPSRO%obs(serialNum:serialNum+systemParameters%varListSize_gpsro-1)%rfict = rfict
        GPSRO%obs(serialNum:serialNum+systemParameters%varListSize_gpsro-1)%z = GPSRO%obs(serialNum)%z
        GPSRO%obs(serialNum:serialNum+systemParameters%varListSize_gpsro-1)%varName = systemParameters%varList_gpsro(1:systemParameters%varListSize_gpsro)
        serialNum = serialNum + (systemParameters%varListSize_gpsro-1)
    enddo
enddo


do iObsVar = 1 , systemParameters%varListSize_gpsro
    where ( adjustl(GPSRO%obs(:)%varName) .eq. adjustl(systemParameters%varList_gpsro(iObsVar)) )
        GPSRO%obs%available = systemParameters%use_varList_gpsro(iObsVar)
    end where
enddo


close(fileID)
!================================================
return
stop
end subroutine getGPSRO
