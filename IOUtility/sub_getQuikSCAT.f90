
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getQuikSCAT(quikscat,systemParameters)

use derivedType

implicit none
type(obsParent),intent(out) :: quikscat
type(systemParameter),intent(in) :: systemParameters

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2

integer ioStatus,fileID
integer i,io,iObsVar  ! loop counter

integer,external :: availableFileID
!================================================

quikscat%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: QuikSCAT'
open(fileID,file='input/obs_input.quikscat',status='old')
do 
    read(fileID,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum
    if ( ioStatus .lt. 0 ) then
        exit
    endif
    quikscat%obsNum = quikscat%obsNum + systemParameters%varListSize_quikscat*levNum
    do i = 1,levNum
        read(fileID,*)  ! just read through lines which are not headers.
    enddo
enddo

! Then, allocate all necessary components.
allocate( quikscat%obs(quikscat%obsNum) )

! Initilize the components just allocated.
quikscat%obs(:)%lon   = 0.0d0
quikscat%obs(:)%lat   = 0.0d0
quikscat%obs(:)%z     = 0.0d0
quikscat%obs(:)%value = 0.0d0
quikscat%obs(:)%error = 0.0d0
quikscat%obs(:)%instrument = 'QuikSCAT'
quikscat%obs(:)%zName      = 'P'
quikscat%obs(:)%varName    = ''
quikscat%obs(:)%available  = .true.
do io = 1 , quikscat%obsNum
    allocate( quikscat%obs(io)%insideHorizontalDomain(systemParameters%max_domain) )
    quikscat%obs(io)%insideHorizontalDomain(:) = .false.
enddo

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) quikscat%obs(serialNum)%z,quikscat%obs(serialNum:serialNum+systemParameters%varListSize_quikscat-1)%value

        quikscat%obs(serialNum:serialNum+systemParameters%varListSize_quikscat-1)%lon = lon
        quikscat%obs(serialNum:serialNum+systemParameters%varListSize_quikscat-1)%lat = lat
        quikscat%obs(serialNum:serialNum+systemParameters%varListSize_quikscat-1)%z   = quikscat%obs(serialNum)%z
        quikscat%obs(serialNum:serialNum+systemParameters%varListSize_quikscat-1)%varName = systemParameters%varList_quikscat(1:systemParameters%varListSize_quikscat)
        serialNum = serialNum + (systemParameters%varListSize_quikscat-1)
    enddo
enddo


quikscat%obs(:)%z = 100.d0 * quikscat%obs(:)%z  !  convert hPa to Pa

do iObsVar = 1 , systemParameters%varListSize_quikscat
    where ( adjustl(quikscat%obs(:)%varName) .eq. adjustl(systemParameters%varList_quikscat(iObsVar)) )
        quikscat%obs%available = systemParameters%use_varList_quikscat(iObsVar)
    end where
enddo


close(fileID)
!================================================
return
stop
end subroutine getQuikSCAT
