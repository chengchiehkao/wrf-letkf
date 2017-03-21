
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getASCAT(ascat,varList,varListSize,use_varList)

use derivedType

implicit none
type(obsParent),intent(out) :: ascat
integer,intent(in)                                   :: varListSize
character(len=10),dimension(varListSize),intent(in)  :: varList
logical,dimension(varListSize),intent(in)            :: use_varList

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2

integer ioStatus,fileID
integer i,iObsVar  ! loop counter

integer,external :: availableFileID
!================================================

ascat%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: ASCAT'
open(fileID,file='input/obs_input.ascat',status='old')
do 
    read(fileID,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum
    if ( ioStatus .lt. 0 ) then
        exit
    endif
    ascat%obsNum = ascat%obsNum + varListSize*levNum
    do i = 1,levNum
        read(fileID,*)  ! just read through lines which are not headers.
    enddo
enddo

! Then, allocate all necessary components.
allocate( ascat%obs(ascat%obsNum) )

! Initilize the components just allocated.
ascat%obs(:)%lon   = 0.0d0
ascat%obs(:)%lat   = 0.0d0
ascat%obs(:)%z     = 0.0d0
ascat%obs(:)%value = 0.0d0
ascat%obs(:)%error = 0.0d0
ascat%obs(:)%instrument = 'ASCAT'
ascat%obs(:)%zName      = 'P'
ascat%obs(:)%varName    = ''
ascat%obs(:)%available  = .true.

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) ascat%obs(serialNum)%z,ascat%obs(serialNum:serialNum+varListSize-1)%value

        ascat%obs(serialNum:serialNum+varListSize-1)%lon = lon
        ascat%obs(serialNum:serialNum+varListSize-1)%lat = lat
        ascat%obs(serialNum:serialNum+varListSize-1)%z   = ascat%obs(serialNum)%z
        ascat%obs(serialNum:serialNum+varListSize-1)%varName = varList(1:varListSize)
        serialNum = serialNum + (varListSize-1)
    enddo
enddo


ascat%obs(:)%z = 100.d0 * ascat%obs(:)%z  !  convert hPa to Pa

do iObsVar = 1 , varListSize
    where ( adjustl(ascat%obs(:)%varName) .eq. adjustl(varList(iObsVar)) )
        ascat%obs%available = use_varList(iObsVar)
    end where
enddo


close(fileID)
!================================================
return
stop
end subroutine getASCAT
