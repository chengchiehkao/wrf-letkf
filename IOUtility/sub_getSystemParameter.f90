
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getSystemParameter(systemParameters)

use derivedType

implicit none

type(systemParameter),intent(out) :: systemParameters

integer :: ensembleSize
logical :: use_sound , use_airep , use_synop , use_amv , use_gpsro , use_airs
integer,parameter :: defaultVarListSize=10
character(len=10),dimension(defaultVarListSize) :: varList_sound , varList_airep , varList_synop , varList_amv , varList_gpsro , varList_airs
logical,dimension(defaultVarListSize) :: use_varList_sound , use_varList_airep , use_varList_synop , use_varList_amv , use_varList_gpsro , use_varList_airs
integer :: varListSize_sound , varListSize_airep , varListSize_synop , varListSize_amv , varListSize_gpsro , varListSize_airs
real(kind=8) :: rd   , rc    ! Decorrelated & highly correlated distance on horizontal space.
real(kind=8) :: rd_z , rc_z  ! Decorrelated & highly correlated distance on vertical space.
real(kind=8) :: inflationFactor

namelist /sizeOfEnsemble/ ensembleSize
namelist /use_observation/ use_sound,use_airep,use_synop,use_amv,use_gpsro,use_airs
namelist /varList_observation/ varList_sound,varList_airep,varList_synop,varList_amv,varList_gpsro,varList_airs
namelist /use_varList_observation/ use_varList_sound,use_varList_airep,use_varList_synop,use_varList_amv,use_varList_gpsro,use_varList_airs
namelist /correlativeDistance/ rd,rc,rd_z,rc_z
namelist /inflation/ inflationFactor

integer :: fileID

integer,external :: availableFileID
!================================================

! Set default value
use_sound = .false.
use_airep = .false.
use_synop = .false.
use_amv   = .false.
use_gpsro = .false.
use_airs  = .false.
use_varList_sound(:) = .true.
use_varList_airep(:) = .true.
use_varList_synop(:) = .true.
use_varList_amv(:)   = .true.
use_varList_gpsro(:) = .true.
use_varList_airs(:)  = .true.
inflationFactor = 1.d0
varList_sound(:) = repeat(' ',len(varList_sound(1)))
varList_airep(:) = repeat(' ',len(varList_airep(1)))
varList_synop(:) = repeat(' ',len(varList_synop(1)))
varList_amv(:)   = repeat(' ',len(varList_amv(1)))
varList_gpsro(:) = repeat(' ',len(varList_gpsro(1)))
varList_airs(:)  = repeat(' ',len(varList_airs(1)))
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
varListSize_airep = count( len_trim(varList_airep) .gt. 0 )
varListSize_synop = count( len_trim(varList_synop) .gt. 0 )
varListSize_amv   = count( len_trim(varList_amv)   .gt. 0 )
varListSize_gpsro = count( len_trim(varList_gpsro) .gt. 0 )
varListSize_airs  = count( len_trim(varList_airs)  .gt. 0 )

if ( use_sound .and. varListSize_sound.eq.0 ) then
    print*,'Variable list of SOUNDING is empty, SOUNDING will be set as disabled.'
    use_sound = .false.
endif
if ( use_airep .and. varListSize_airep.eq.0 ) then
    print*,'Variable list of AIREP is empty, AIREP will be set as disabled.'
    use_airep = .false.
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
if ( use_airs .and. varListSize_airs.eq.0 ) then
    print*,'Variable list of AIRS is empty, AIRS will be set as disabled.'
    use_airs = .false.
endif


if ( rd .le. 0.d0 ) then
    print*,'Rd shall > 0., program stopped.'
    stop
endif
if ( rc .le. 0.d0 ) then
    print*,'Rc shall > 0., program stopped.'
    stop
endif
if ( rd_z .le. 0.d0 ) then
    print*,'Rd_z shall > 0., program stopped.'
    stop
endif
if ( rc_z .le. 0.d0 ) then
    print*,'Rc_z shall > 0., program stopped.'
    stop
endif

if ( inflationFactor .le. 0.d0 ) then
    print*,'inflation factor shall > 0., will automatically set to be 1.'
    inflationFactor = 1.d0
endif


if ( use_sound )  allocate( systemParameters % varList_sound(varListSize_sound) )
if ( use_airep )  allocate( systemParameters % varList_airep(varListSize_airep) )
if ( use_synop )  allocate( systemParameters % varList_synop(varListSize_synop) )
if ( use_amv   )  allocate( systemParameters % varList_amv(varListSize_amv) )
if ( use_gpsro )  allocate( systemParameters % varList_gpsro(varListSize_gpsro) )
if ( use_airs )   allocate( systemParameters % varList_airs(varListSize_airs) )

if ( use_sound )  allocate( systemParameters % use_varList_sound(varListSize_sound) )
if ( use_airep )  allocate( systemParameters % use_varList_airep(varListSize_airep) )
if ( use_synop )  allocate( systemParameters % use_varList_synop(varListSize_synop) )
if ( use_amv   )  allocate( systemParameters % use_varList_amv(varListSize_amv) )
if ( use_gpsro )  allocate( systemParameters % use_varList_gpsro(varListSize_gpsro) )
if ( use_airs  )  allocate( systemParameters % use_varList_airs(varListSize_airs) )

if ( use_sound )  systemParameters % varList_sound(:) = repeat(' ',len(varList_sound(1)))
if ( use_airep )  systemParameters % varList_airep(:) = repeat(' ',len(varList_airep(1)))
if ( use_synop )  systemParameters % varList_synop(:) = repeat(' ',len(varList_synop(1)))
if ( use_amv   )  systemParameters % varList_amv(:)   = repeat(' ',len(varList_amv(1)))
if ( use_gpsro )  systemParameters % varList_gpsro(:) = repeat(' ',len(varList_gpsro(1)))
if ( use_airs )   systemParameters % varList_airs(:)  = repeat(' ',len(varList_airs(1)))


systemParameters % ensembleSize = ensembleSize
systemParameters % use_sound = use_sound
systemParameters % use_airep = use_airep
systemParameters % use_synop = use_synop
systemParameters % use_amv   = use_amv
systemParameters % use_gpsro = use_gpsro
systemParameters % use_airs  = use_airs
systemParameters % varListSize_sound = varListSize_sound
systemParameters % varListSize_airep = varListSize_airep
systemParameters % varListSize_synop = varListSize_synop
systemParameters % varListSize_amv   = varListSize_amv
systemParameters % varListSize_gpsro = varListSize_gpsro
systemParameters % varListSize_airs  = varListSize_airs
if ( use_sound )  systemParameters % varList_sound(:) = varList_sound(1:varListSize_sound)
if ( use_airep )  systemParameters % varList_airep(:) = varList_airep(1:varListSize_airep)
if ( use_synop )  systemParameters % varList_synop(:) = varList_synop(1:varListSize_synop)
if ( use_amv   )  systemParameters % varList_amv(:)   = varList_amv(1:varListSize_amv)
if ( use_gpsro )  systemParameters % varList_gpsro(:) = varList_gpsro(1:varListSize_gpsro)
if ( use_airs  )  systemParameters % varList_airs(:)  = varList_airs(1:varListSize_airs)

if ( use_sound )  systemParameters % use_varList_sound(:) = use_varList_sound(1:varListSize_sound)
if ( use_airep )  systemParameters % use_varList_airep(:) = use_varList_airep(1:varListSize_airep)
if ( use_synop )  systemParameters % use_varList_synop(:) = use_varList_synop(1:varListSize_synop)
if ( use_amv   )  systemParameters % use_varList_amv(:)   = use_varList_amv(1:varListSize_amv)
if ( use_gpsro )  systemParameters % use_varList_gpsro(:) = use_varList_gpsro(1:varListSize_gpsro)
if ( use_airs  )  systemParameters % use_varList_airs(:)  = use_varList_airs(1:varListSize_airs)
systemParameters % rd = rd
systemParameters % rc = rc
systemParameters % rd_z = rd_z
systemParameters % rc_z = rc_z
systemParameters % inflationFactor = inflationFactor


!================================================
return
stop
end subroutine getSystemParameter
