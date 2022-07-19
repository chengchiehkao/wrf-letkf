
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getIASI(iasi,systemParameters)

use derivedType

implicit none
type(obsParent),intent(out) :: iasi
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

iasi%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: IASI'
open(fileID,file='input/obs_input.IASI',status='old')
do 
    header = ''
    read(fileID,'(a100)',iostat=ioStatus) header
    if ( ioStatus .lt. 0 ) then
        exit
    endif

    read(header,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum

    iasi%obsNum = iasi%obsNum + systemParameters%varListSize_iasi*levNum
    do i = 1,levNum
        read(fileID,*)  ! just read through lines which are not headers.
    enddo

enddo

! Then, allocate all necessary components.
allocate( iasi%obs(iasi%obsNum) )

! Initilize the components just allocated.
iasi%obs(:)%lon   = 0.0d0
iasi%obs(:)%lat   = 0.0d0
iasi%obs(:)%z     = 0.0d0
iasi%obs(:)%value = 0.0d0
iasi%obs(:)%error = 0.0d0
iasi%obs(:)%instrument = 'IASI'
iasi%obs(:)%zName      = 'P'
iasi%obs(:)%varName    = ''
iasi%obs(:)%available  = .true.
do io = 1 , iasi%obsNum
    allocate( iasi%obs(io)%insideHorizontalDomain(systemParameters%max_domain) )
    iasi%obs(io)%insideHorizontalDomain(:) = .false.
enddo

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
        read(fileID,*) iasi%obs(serialNum)%z,iasi%obs(serialNum:serialNum+systemParameters%varListSize_iasi-1)%value

        iasi%obs(serialNum:serialNum+systemParameters%varListSize_iasi-1)%lon = lon
        iasi%obs(serialNum:serialNum+systemParameters%varListSize_iasi-1)%lat = lat
        iasi%obs(serialNum:serialNum+systemParameters%varListSize_iasi-1)%z   = iasi%obs(serialNum)%z
        iasi%obs(serialNum:serialNum+systemParameters%varListSize_iasi-1)%varName = systemParameters%varList_iasi(1:systemParameters%varListSize_iasi)
        serialNum = serialNum + (systemParameters%varListSize_iasi-1)
    enddo

enddo


! Convert temperature in Celsius to Kelvin
do io = 1 , iasi%obsNum
    if ( trim(adjustl(iasi%obs(io)%varName)).eq.'T' .and. &
         dabs(iasi%obs(io)%value - invalidValue) .gt. dabs(invalidValue*epsilon(1.d0)) ) then
        iasi%obs(io)%value = iasi%obs(io)%value + 273.15d0
    endif
enddo


iasi%obs(:)%z = 100.d0 * iasi%obs(:)%z  !  convert hPa to Pa

do iObsVar = 1 , systemParameters%varListSize_iasi
    where ( adjustl(iasi%obs(:)%varName) .eq. adjustl(systemParameters%varList_iasi(iObsVar)) )
        iasi%obs%available = systemParameters%use_varList_iasi(iObsVar)
    end where
enddo


close(fileID)
!================================================
return
stop
end subroutine getIASI
