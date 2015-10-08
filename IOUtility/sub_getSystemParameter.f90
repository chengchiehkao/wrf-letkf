
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getSystemParameter(systemParameters)

use derivedType

implicit none

type(systemParameter),intent(out) :: systemParameters

integer :: ensembleSize
logical :: use_sound , use_synop , use_amv , use_gpsro
integer,parameter :: defaultVarListSize=10
character(len=10),dimension(defaultVarListSize) :: varList_sound , varList_synop , varList_amv , varList_gpsro
logical,dimension(defaultVarListSize) :: use_varList_sound , use_varList_synop , use_varList_amv , use_varList_gpsro
integer :: varListSize_sound , varListSize_synop , varListSize_amv , varListSize_gpsro
real(kind=8) :: rd , rc  ! Decorrelated & highly correlated distance on horizontal space.
real(kind=8) :: rv       ! Decorrelated distance on vertical space.teger ioStatus,fileID
real(kind=8) :: inflationFactor

namelist /sizeOfEnsemble/ ensembleSize
namelist /use_observation/ use_sound,use_synop,use_amv,use_gpsro
namelist /varList_observation/ varList_sound,varList_synop,varList_amv,varList_gpsro
namelist /use_varList_observation/ use_varList_sound,use_varList_synop,use_varList_amv,use_varList_gpsro
namelist /correlativeDistance/ rd,rc,rv
namelist /inflation/ inflationFactor

integer :: fileID

integer,external :: availableFileID
!================================================

! Set default value
use_sound = .false.
use_synop = .false.
use_amv   = .false.
use_gpsro = .false.
use_varList_sound(:) = .true.
use_varList_synop(:) = .true.
use_varList_amv(:)   = .true.
use_varList_gpsro(:) = .true.
inflationFactor = 1.d0
varList_sound(:) = repeat(' ',len(varList_sound(1)))
varList_synop(:) = repeat(' ',len(varList_synop(1)))
varList_amv(:)   = repeat(' ',len(varList_amv(1)))
varList_gpsro(:) = repeat(' ',len(varList_gpsro(1)))
! End of assignment


fileID = availableFileID()
open(fileID,file='input/systemParameter.nml',status='old')
read(fileID,nml = sizeOfEnsemble)
read(fileID,nml = correlativeDistance)
read(fileID,nml = inflation)
read(fileID,nml = use_observation)
read(fileID,nml = varList_observation)
read(fileID,nml = use_varList_observation)
close(fileID)

if ( ensembleSize .lt. 2 ) then
    print*,'Ensemble size shall >= 2, program stopped.'
    stop
endif

! Shall beware of non-consecutive variable list.
varListSize_sound = count( len_trim(varList_sound) .gt. 0 )
varListSize_synop = count( len_trim(varList_synop) .gt. 0 )
varListSize_amv   = count( len_trim(varList_amv)   .gt. 0 )
varListSize_gpsro = count( len_trim(varList_gpsro) .gt. 0 )

if ( use_sound .and. varListSize_sound.eq.0 ) then
    print*,'Variable list of SOUNDING is empty, SOUNDING will be set as disabled.'
    use_sound = .false.
endif
if ( use_synop .and. varListSize_synop.eq.0 ) then
    print*,'Variable list of SYNOP is empty, SYNOP will be set as disabled.'
    use_synop = .false.
endif
if ( use_amv .and. varListSize_amv.eq.0 ) then
    print*,'Variable list of AMV is empty, AMV will be set as disabled.'
    use_amv = .false.
endif
if ( use_gpsro .and. varListSize_gpsro.eq.0 ) then
    print*,'Variable list of GPSRO is empty, GPSRO will be set as disabled.'
    use_gpsro = .false.
endif


if ( rd .le. 0.d0 ) then
    print*,'Rd shall > 0., program stopped.'
    stop
endif
if ( rc .le. 0.d0 ) then
    print*,'Rc shall > 0., program stopped.'
    stop
endif
if ( rv .le. 0.d0 ) then
    print*,'Rv shall > 0., program stopped.'
    stop
endif

if ( inflationFactor .le. 0.d0 ) then
    print*,'inflation factor shall > 0., will automatically set to be 1.'
    inflationFactor = 1.d0
endif


if ( use_sound )  allocate( systemParameters % varList_sound(varListSize_sound) )
if ( use_synop )  allocate( systemParameters % varList_synop(varListSize_synop) )
if ( use_amv   )  allocate( systemParameters % varList_amv(varListSize_amv) )
if ( use_gpsro )  allocate( systemParameters % varList_gpsro(varListSize_gpsro) )

if ( use_sound )  allocate( systemParameters % use_varList_sound(varListSize_sound) )
if ( use_synop )  allocate( systemParameters % use_varList_synop(varListSize_synop) )
if ( use_amv   )  allocate( systemParameters % use_varList_amv(varListSize_amv) )
if ( use_gpsro )  allocate( systemParameters % use_varList_gpsro(varListSize_gpsro) )

if ( use_sound )  systemParameters % varList_sound(:) = repeat(' ',len(varList_sound(1)))
if ( use_synop )  systemParameters % varList_synop(:) = repeat(' ',len(varList_synop(1)))
if ( use_amv   )  systemParameters % varList_amv(:)   = repeat(' ',len(varList_amv(1)))
if ( use_gpsro )  systemParameters % varList_gpsro(:) = repeat(' ',len(varList_gpsro(1)))


systemParameters % ensembleSize = ensembleSize
systemParameters % use_sound = use_sound
systemParameters % use_synop = use_synop
systemParameters % use_amv   = use_amv
systemParameters % use_gpsro = use_gpsro
systemParameters % varListSize_sound = varListSize_sound
systemParameters % varListSize_synop = varListSize_synop
systemParameters % varListSize_amv   = varListSize_amv
systemParameters % varListSize_gpsro = varListSize_gpsro
if ( use_sound )  systemParameters % varList_sound(:) = varList_sound(1:varListSize_sound)
if ( use_synop )  systemParameters % varList_synop(:) = varList_synop(1:varListSize_synop)
if ( use_amv   )  systemParameters % varList_amv(:)   = varList_amv(1:varListSize_amv)
if ( use_gpsro )  systemParameters % varList_gpsro(:) = varList_gpsro(1:varListSize_gpsro)

if ( use_sound )  systemParameters % use_varList_sound(:) = use_varList_sound(1:varListSize_sound)
if ( use_synop )  systemParameters % use_varList_synop(:) = use_varList_synop(1:varListSize_synop)
if ( use_amv   )  systemParameters % use_varList_amv(:)   = use_varList_amv(1:varListSize_amv)
if ( use_gpsro )  systemParameters % use_varList_gpsro(:) = use_varList_gpsro(1:varListSize_gpsro)
systemParameters % rd = rd
systemParameters % rc = rc
systemParameters % rv = rv
systemParameters % inflationFactor = inflationFactor


!================================================
return
stop
end subroutine getSystemParameter
