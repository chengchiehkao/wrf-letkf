
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getAMV(amv,systemParameters)

use derivedType

implicit none
type(obsParent),intent(out) :: amv
type(systemParameter),intent(in) :: systemParameters

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2

integer ioStatus,fileID
integer i,io,iObsVar  ! loop counter

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
    amv%obsNum = amv%obsNum + systemParameters%varListSize_amv*levNum
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
do io = 1 , amv%obsNum
    allocate( amv%obs(io)%insideHorizontalDomain(systemParameters%max_domain) )
    amv%obs(io)%insideHorizontalDomain(:) = .false.
enddo

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) amv%obs(serialNum)%z,amv%obs(serialNum:serialNum+systemParameters%varListSize_amv-1)%value

        amv%obs(serialNum:serialNum+systemParameters%varListSize_amv-1)%lon = lon
        amv%obs(serialNum:serialNum+systemParameters%varListSize_amv-1)%lat = lat
        amv%obs(serialNum:serialNum+systemParameters%varListSize_amv-1)%z   = amv%obs(serialNum)%z
        amv%obs(serialNum:serialNum+systemParameters%varListSize_amv-1)%varName = systemParameters%varList_amv(1:systemParameters%varListSize_amv)
        serialNum = serialNum + (systemParameters%varListSize_amv-1)
    enddo
enddo


amv%obs(:)%z = 100.d0 * amv%obs(:)%z  !  convert hPa to Pa

do iObsVar = 1 , systemParameters%varListSize_amv
    where ( adjustl(amv%obs(:)%varName) .eq. adjustl(systemParameters%varList_amv(iObsVar)) )
        amv%obs%available = systemParameters%use_varList_amv(iObsVar)
    end where
enddo


close(fileID)
!================================================
return
stop
end subroutine getAMV
