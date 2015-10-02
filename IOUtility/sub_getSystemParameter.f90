
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
integer :: varListSize_sound , varListSize_synop , varListSize_amv , varListSize_gpsro
real(kind=8) :: rd , rc  ! Decorrelated & highly correlated distance on horizontal space.
real(kind=8) :: rv       ! Decorrelated distance on vertical space.teger ioStatus,fileID

namelist /sizeOfEnsemble/ ensembleSize
namelist /use_observation/ use_sound,use_synop,use_amv,use_gpsro
namelist /varList_observation/ varList_sound,varList_synop,varList_amv,varList_gpsro
namelist /correlativeDistance/ rd,rc,rv

integer :: fileID

integer,external :: availableFileID
!================================================

fileID = availableFileID()
open(fileID,file='systemParameter.nml',status='old')
read(fileID,nml = sizeOfEnsemble)
read(fileID,nml = use_observation)
read(fileID,nml = varList_observation)
read(fileID,nml = correlativeDistance)
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


if ( use_sound )  allocate( systemParameters % varList_sound(varListSize_sound) )
if ( use_synop )  allocate( systemParameters % varList_synop(varListSize_synop) )
if ( use_amv   )  allocate( systemParameters % varList_amv(varListSize_amv) )
if ( use_gpsro )  allocate( systemParameters % varList_gpsro(varListSize_gpsro) )

if ( use_sound )  systemParameters % varList_sound(:) = repeat(' ',len(varList_sound))
if ( use_synop )  systemParameters % varList_synop(:) = repeat(' ',len(varList_synop))
if ( use_amv   )  systemParameters % varList_amv(:)   = repeat(' ',len(varList_amv))
if ( use_gpsro )  systemParameters % varList_gpsro(:) = repeat(' ',len(varList_gpsro))


systemParameters % ensembleSize =  ensembleSize
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
systemParameters % rd = rd
systemParameters % rc = rc
systemParameters % rv = rv


!================================================
return
stop
end subroutine getSystemParameter
