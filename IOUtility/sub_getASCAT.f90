
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getASCAT(ascat,systemParameters)

use derivedType

implicit none
type(obsParent),intent(out) :: ascat
type(systemParameter),intent(in) :: systemParameters

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2

integer ioStatus,fileID
integer i,io,iObsVar  ! loop counter

integer,external :: availableFileID
!================================================

ascat%obsNum = 0  ; print*,'entering getASCAT'

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: ASCAT'
open(fileID,file='input/obs_input.ascat',status='old')
do 
    read(fileID,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum
    if ( ioStatus .lt. 0 ) then
        exit
    endif
    ascat%obsNum = ascat%obsNum + systemParameters%varListSize_ascat*levNum
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
do io = 1 , ascat%obsNum
    allocate( ascat%obs(io)%insideHorizontalDomain(systemParameters%max_domain) )
    ascat%obs(io)%insideHorizontalDomain(:) = .false.
enddo

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) ascat%obs(serialNum)%z,ascat%obs(serialNum:serialNum+systemParameters%varListSize_ascat-1)%value

        ascat%obs(serialNum:serialNum+systemParameters%varListSize_ascat-1)%lon = lon
        ascat%obs(serialNum:serialNum+systemParameters%varListSize_ascat-1)%lat = lat
        ascat%obs(serialNum:serialNum+systemParameters%varListSize_ascat-1)%z   = ascat%obs(serialNum)%z
        ascat%obs(serialNum:serialNum+systemParameters%varListSize_ascat-1)%varName = systemParameters%varList_ascat(1:systemParameters%varListSize_ascat)
        serialNum = serialNum + (systemParameters%varListSize_ascat-1)
    enddo
enddo


ascat%obs(:)%z = 100.d0 * ascat%obs(:)%z  !  convert hPa to Pa

do iObsVar = 1 , systemParameters%varListSize_ascat
    where ( adjustl(ascat%obs(:)%varName) .eq. adjustl(systemParameters%varList_ascat(iObsVar)) )
        ascat%obs%available = systemParameters%use_varList_ascat(iObsVar)
    end where
enddo


close(fileID)  ; print*,'leaving getASCAT'
!================================================
return
stop
end subroutine getASCAT
