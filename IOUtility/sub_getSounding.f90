
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getSounding(sounding,varList,varListSize)

use derivedType

implicit none
type(obsParent),intent(out) :: sounding
integer,intent(in)                                   :: varListSize
character(len=10),dimension(varListSize),intent(in)  :: varList

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2
character(len=100) header
real(kind=8),parameter :: invalidValue = -888888.d0

integer ioStatus,fileID
integer i,io  ! loop counter

integer,external :: availableFileID
!================================================

sounding%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: sounding'
open(fileID,file='input/obs_input.sounding',status='old')
do 
    header = ''
    read(fileID,'(a100)',iostat=ioStatus) header
    if ( ioStatus .lt. 0 ) then
        exit
    endif

    read(header,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum

    if ( index(header,'AIREP') .eq. 0 ) then  ! means it's not AIREP
        sounding%obsNum = sounding%obsNum + varListSize*levNum +1  ! +1 for PSFC
        do i = 1,levNum+1  ! +1 for PSFC
            read(fileID,*)  ! just read through lines which are not headers.
        enddo
    elseif ( index(header,'AIREP') .ne. 0 ) then  ! means it's AIREP
        sounding%obsNum = sounding%obsNum + varListSize*levNum
        do i = 1,levNum
            read(fileID,*)  ! just read through lines which are not headers.
        enddo
    endif
enddo

! Then, allocate all necessary components.
allocate( sounding%obs(sounding%obsNum) )

! Initilize the components just allocated.
sounding%obs(:)%lon   = 0.0d0
sounding%obs(:)%lat   = 0.0d0
sounding%obs(:)%z     = 0.0d0
sounding%obs(:)%value = 0.0d0
sounding%obs(:)%error = 0.0d0
sounding%obs(:)%type       = 'SOUNDING  '
sounding%obs(:)%instrument = 'SOUNDING  '
sounding%obs(:)%zName      = ''
sounding%obs(:)%varName    = ''
sounding%obs(:)%available  = .true.

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    header = ''
    read(fileID,'(a100)',iostat=ioStatus) header
    if ( ioStatus .lt. 0 ) exit
    read(header,*) lon,lat,levNum

    if ( index(header,'AIREP') .eq. 0 ) then  ! means it's not AIREP
        serialNum = serialNum + 1
        sounding%obs(serialNum)%lon = lon
        sounding%obs(serialNum)%lat = lat
        sounding%obs(serialNum)%zName = 'SURFACE   '
        sounding%obs(serialNum)%varName = 'PSFC      '
        read(fileID,*) sounding%obs(serialNum)%value
        do i = 1,levNum
            serialNum = serialNum + 1
            read(fileID,*) sounding%obs(serialNum)%z,sounding%obs(serialNum:serialNum+varListSize-1)%value

            sounding%obs(serialNum:serialNum+varListSize-1)%lon = lon
            sounding%obs(serialNum:serialNum+varListSize-1)%lat = lat
            sounding%obs(serialNum:serialNum+varListSize-1)%z   = sounding%obs(serialNum)%z
            sounding%obs(serialNum:serialNum+varListSize-1)%zName = 'P         '
            sounding%obs(serialNum:serialNum+varListSize-1)%varName = varList(1:varListSize)
            serialNum = serialNum + (varListSize-1)
        enddo
    elseif ( index(header,'AIREP') .ne. 0 ) then  ! means it's AIREP
        do i = 1,levNum
            serialNum = serialNum + 1
            read(fileID,*) sounding%obs(serialNum)%z,sounding%obs(serialNum:serialNum+varListSize-1)%value

            sounding%obs(serialNum:serialNum+varListSize-1)%instrument = 'AIREP     '
            sounding%obs(serialNum:serialNum+varListSize-1)%lon = lon
            sounding%obs(serialNum:serialNum+varListSize-1)%lat = lat
            sounding%obs(serialNum:serialNum+varListSize-1)%z   = sounding%obs(serialNum)%z
            sounding%obs(serialNum:serialNum+varListSize-1)%zName = 'P         '
            sounding%obs(serialNum:serialNum+varListSize-1)%varName = varList(1:varListSize)
            serialNum = serialNum + (varListSize-1)
        enddo
    endif
enddo


! Convert temperature in Celsius to Kelvin
!where ( trim(adjustl(sounding%obs(:)%varName)).eq.'T' .and. &
!        dabs(sounding%obs(:)%var - invalidValue) .gt. dabs(invalidValue*epsilon(1.d0)) )
!
!    sounding%obs(:)%var = sounding%obs(:)%var + 273.d0  ! not 273.15 but 273. because of the preprocess of sounding obs
!end where

do io = 1 , sounding%obsNum
    if ( trim(adjustl(sounding%obs(io)%varName)).eq.'T' .and. &
         dabs(sounding%obs(io)%value - invalidValue) .gt. dabs(invalidValue*epsilon(1.d0)) ) then
        sounding%obs(io)%value = sounding%obs(io)%value + 273.15d0
    endif
enddo

do io = 1 , sounding%obsNum
    if ( trim(adjustl(sounding%obs(io)%zName)) .ne. 'SURFACE' ) then
        sounding%obs(io)%z = 100.d0 * sounding%obs(io)%z  !  convert hPa to Pa
    endif
enddo


close(fileID)
!================================================
return
stop
end subroutine getSounding
