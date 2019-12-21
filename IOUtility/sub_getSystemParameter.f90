
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getSystemParameter(systemParameters)

use derivedType

implicit none

type(systemParameter),intent(out) :: systemParameters

integer :: ensembleSize
logical :: use_sound , use_airep , use_synop , use_amv , use_gpsro , use_airs , use_quikscat , use_ascat , use_iasi , use_oscat , use_windsat , use_cygnss
integer,parameter :: defaultVarListSize=10
character(len=10),dimension(defaultVarListSize) :: varList_sound , varList_airep , varList_synop , varList_amv , varList_gpsro , varList_airs , varList_quikscat , varList_ascat , varList_iasi , varList_oscat , varList_windsat , varList_cygnss
logical,dimension(defaultVarListSize) :: use_varList_sound , use_varList_airep , use_varList_synop , use_varList_amv , use_varList_gpsro , &
                                         use_varList_airs , use_varList_quikscat , use_varList_ascat , use_varList_iasi , use_varList_oscat , &
                                         use_varList_windsat , use_varList_cygnss
integer :: varListSize_sound , varListSize_airep , varListSize_synop , varListSize_amv , varListSize_gpsro , varListSize_airs , varListSize_quikscat , varListSize_ascat , varListSize_iasi , varListSize_oscat , varListSize_windsat , varListSize_cygnss
real(kind=8) :: rd   , rc    ! Decorrelated & highly correlated distance on horizontal space.
real(kind=8) :: rd_z , rc_z  ! Decorrelated & highly correlated distance on vertical space.
real(kind=8) :: inflationFactor
integer :: boundaryWidth

namelist /sizeOfEnsemble/ ensembleSize
namelist /domain/ boundaryWidth
namelist /use_observation/ use_sound,use_airep,use_synop,use_amv,use_gpsro,use_airs,use_quikscat,use_ascat,use_iasi,use_oscat,use_windsat,use_cygnss
namelist /varList_observation/ varList_sound,varList_airep,varList_synop,varList_amv,varList_gpsro,varList_airs,varList_quikscat,varList_ascat,varList_iasi,varList_oscat,varList_windsat,varList_cygnss
namelist /use_varList_observation/ use_varList_sound,use_varList_airep,use_varList_synop,use_varList_amv,use_varList_gpsro,use_varList_airs,use_varList_quikscat,use_varList_ascat,use_varList_iasi,use_varList_oscat,use_varList_windsat,use_varList_cygnss
namelist /correlativeDistance/ rd,rc,rd_z,rc_z
namelist /inflation/ inflationFactor

integer :: fileID

integer,external :: availableFileID
!================================================

! Set default value
boundaryWidth = 5
use_sound    = .false.
use_airep    = .false.
use_synop    = .false.
use_amv      = .false.
use_gpsro    = .false.
use_airs     = .false.
use_quikscat = .false.
use_ascat    = .false.
use_iasi     = .false.
use_oscat    = .false.
use_windsat  = .false.
use_cygnss   = .false.
use_varList_sound(:)    = .true.
use_varList_airep(:)    = .true.
use_varList_synop(:)    = .true.
use_varList_amv(:)      = .true.
use_varList_gpsro(:)    = .true.
use_varList_airs(:)     = .true.
use_varList_quikscat(:) = .true.
use_varList_ascat(:)    = .true.
use_varList_iasi(:)     = .true.
use_varList_oscat(:)    = .true.
use_varList_windsat(:)  = .true.
use_varList_cygnss(:)   = .true.
inflationFactor = 1.d0
varList_sound(:)    = repeat(' ',len(varList_sound(1)))
varList_airep(:)    = repeat(' ',len(varList_airep(1)))
varList_synop(:)    = repeat(' ',len(varList_synop(1)))
varList_amv(:)      = repeat(' ',len(varList_amv(1)))
varList_gpsro(:)    = repeat(' ',len(varList_gpsro(1)))
varList_airs(:)     = repeat(' ',len(varList_airs(1)))
varList_quikscat(:) = repeat(' ',len(varList_quikscat(1)))
varList_ascat(:)    = repeat(' ',len(varList_ascat(1)))
varList_iasi(:)     = repeat(' ',len(varList_iasi(1)))
varList_oscat(:)    = repeat(' ',len(varList_oscat(1)))
varList_windsat(:)  = repeat(' ',len(varList_windsat(1)))
varList_cygnss(:)   = repeat(' ',len(varList_cygnss(1)))
! End of assignment


fileID = availableFileID()
open(fileID,file='input/systemParameter.nml',status='old')
read(fileID,nml = sizeOfEnsemble)
read(fileID,nml = domain)
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

if ( boundaryWidth .lt. 0 ) then
    print*,'Boundary width shall >= 0, program stopped.'
    stop
endif

! Shall beware of non-consecutive variable list.
varListSize_sound     = count( len_trim(varList_sound) .gt. 0 )
varListSize_airep     = count( len_trim(varList_airep) .gt. 0 )
varListSize_synop     = count( len_trim(varList_synop) .gt. 0 )
varListSize_amv       = count( len_trim(varList_amv)   .gt. 0 )
varListSize_gpsro     = count( len_trim(varList_gpsro) .gt. 0 )
varListSize_airs      = count( len_trim(varList_airs)  .gt. 0 )
varListSize_quikscat  = count( len_trim(varList_quikscat)  .gt. 0 )
varListSize_ascat     = count( len_trim(varList_ascat) .gt. 0 )
varListSize_iasi      = count( len_trim(varList_iasi)  .gt. 0 )
varListSize_oscat     = count( len_trim(varList_oscat) .gt. 0 )
varListSize_windsat   = count( len_trim(varList_windsat) .gt. 0 )
varListSize_cygnss    = count( len_trim(varList_cygnss) .gt. 0 )

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
if ( use_quikscat .and. varListSize_quikscat.eq.0 ) then
    print*,'Variable list of QuikSCAT is empty, QuikSCAT will be set as disabled.'
    use_quikscat = .false.
endif
if ( use_ascat .and. varListSize_ascat.eq.0 ) then
    print*,'Variable list of ASCAT is empty, ASCAT will be set as disabled.'
    use_ascat = .false.
endif
if ( use_iasi .and. varListSize_iasi.eq.0 ) then
    print*,'Variable list of IASI is empty, IASI will be set as disabled.'
    use_iasi = .false.
endif
if ( use_oscat .and. varListSize_oscat.eq.0 ) then
    print*,'Variable list of OSCAT is empty, OSCAT will be set as disabled.'
    use_oscat = .false.
endif
if ( use_windsat .and. varListSize_windsat.eq.0 ) then
    print*,'Variable list of WindSat is empty, WindSat will be set as disabled.'
    use_windsat = .false.
endif
if ( use_cygnss .and. varListSize_cygnss.eq.0 ) then
    print*,'Variable list of CYGNSS is empty, CYGNSS will be set as disabled.'
    use_cygnss = .false.
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


if ( use_sound    )  allocate( systemParameters % varList_sound(varListSize_sound) )
if ( use_airep    )  allocate( systemParameters % varList_airep(varListSize_airep) )
if ( use_synop    )  allocate( systemParameters % varList_synop(varListSize_synop) )
if ( use_amv      )  allocate( systemParameters % varList_amv(varListSize_amv) )
if ( use_gpsro    )  allocate( systemParameters % varList_gpsro(varListSize_gpsro) )
if ( use_airs     )  allocate( systemParameters % varList_airs(varListSize_airs) )
if ( use_quikscat )  allocate( systemParameters % varList_quikscat(varListSize_quikscat) )
if ( use_ascat    )  allocate( systemParameters % varList_ascat(varListSize_ascat) )
if ( use_iasi     )  allocate( systemParameters % varList_iasi(varListSize_iasi) )
if ( use_oscat    )  allocate( systemParameters % varList_oscat(varListSize_oscat) )
if ( use_windsat  )  allocate( systemParameters % varList_windsat(varListSize_windsat) )
if ( use_cygnss   )  allocate( systemParameters % varList_cygnss(varListSize_cygnss) )

if ( use_sound    )  allocate( systemParameters % use_varList_sound(varListSize_sound) )
if ( use_airep    )  allocate( systemParameters % use_varList_airep(varListSize_airep) )
if ( use_synop    )  allocate( systemParameters % use_varList_synop(varListSize_synop) )
if ( use_amv      )  allocate( systemParameters % use_varList_amv(varListSize_amv) )
if ( use_gpsro    )  allocate( systemParameters % use_varList_gpsro(varListSize_gpsro) )
if ( use_airs     )  allocate( systemParameters % use_varList_airs(varListSize_airs) )
if ( use_quikscat )  allocate( systemParameters % use_varList_quikscat(varListSize_quikscat) )
if ( use_ascat    )  allocate( systemParameters % use_varList_ascat(varListSize_ascat) )
if ( use_iasi     )  allocate( systemParameters % use_varList_iasi(varListSize_iasi) )
if ( use_oscat    )  allocate( systemParameters % use_varList_oscat(varListSize_oscat) )
if ( use_windsat  )  allocate( systemParameters % use_varList_windsat(varListSize_windsat) )
if ( use_cygnss   )  allocate( systemParameters % use_varList_cygnss(varListSize_cygnss) )

if ( use_sound    )  systemParameters % varList_sound(:)    = repeat(' ',len(varList_sound(1)))
if ( use_airep    )  systemParameters % varList_airep(:)    = repeat(' ',len(varList_airep(1)))
if ( use_synop    )  systemParameters % varList_synop(:)    = repeat(' ',len(varList_synop(1)))
if ( use_amv      )  systemParameters % varList_amv(:)      = repeat(' ',len(varList_amv(1)))
if ( use_gpsro    )  systemParameters % varList_gpsro(:)    = repeat(' ',len(varList_gpsro(1)))
if ( use_airs     )  systemParameters % varList_airs(:)     = repeat(' ',len(varList_airs(1)))
if ( use_quikscat )  systemParameters % varList_quikscat(:) = repeat(' ',len(varList_quikscat(1)))
if ( use_ascat    )  systemParameters % varList_ascat(:)    = repeat(' ',len(varList_ascat(1)))
if ( use_iasi     )  systemParameters % varList_iasi(:)     = repeat(' ',len(varList_iasi(1)))
if ( use_oscat    )  systemParameters % varList_oscat(:)    = repeat(' ',len(varList_oscat(1)))
if ( use_windsat  )  systemParameters % varList_windsat(:)  = repeat(' ',len(varList_windsat(1)))
if ( use_cygnss   )  systemParameters % varList_cygnss(:)    = repeat(' ',len(varList_cygnss(1)))


systemParameters % ensembleSize = ensembleSize
systemParameters % boundaryWidth = boundaryWidth
systemParameters % use_sound    = use_sound
systemParameters % use_airep    = use_airep
systemParameters % use_synop    = use_synop
systemParameters % use_amv      = use_amv
systemParameters % use_gpsro    = use_gpsro
systemParameters % use_airs     = use_airs
systemParameters % use_quikscat = use_quikscat
systemParameters % use_ascat    = use_ascat
systemParameters % use_iasi     = use_iasi
systemParameters % use_oscat    = use_oscat
systemParameters % use_windsat  = use_windsat
systemParameters % use_cygnss   = use_cygnss
systemParameters % varListSize_sound    = varListSize_sound
systemParameters % varListSize_airep    = varListSize_airep
systemParameters % varListSize_synop    = varListSize_synop
systemParameters % varListSize_amv      = varListSize_amv
systemParameters % varListSize_gpsro    = varListSize_gpsro
systemParameters % varListSize_airs     = varListSize_airs
systemParameters % varListSize_quikscat = varListSize_quikscat
systemParameters % varListSize_ascat    = varListSize_ascat
systemParameters % varListSize_iasi     = varListSize_iasi
systemParameters % varListSize_oscat    = varListSize_oscat
systemParameters % varListSize_windsat  = varListSize_windsat
systemParameters % varListSize_cygnss   = varListSize_cygnss
if ( use_sound    )  systemParameters % varList_sound(:)    = varList_sound(1:varListSize_sound)
if ( use_airep    )  systemParameters % varList_airep(:)    = varList_airep(1:varListSize_airep)
if ( use_synop    )  systemParameters % varList_synop(:)    = varList_synop(1:varListSize_synop)
if ( use_amv      )  systemParameters % varList_amv(:)      = varList_amv(1:varListSize_amv)
if ( use_gpsro    )  systemParameters % varList_gpsro(:)    = varList_gpsro(1:varListSize_gpsro)
if ( use_airs     )  systemParameters % varList_airs(:)     = varList_airs(1:varListSize_airs)
if ( use_quikscat )  systemParameters % varList_quikscat(:) = varList_quikscat(1:varListSize_quikscat)
if ( use_ascat    )  systemParameters % varList_ascat(:)    = varList_ascat(1:varListSize_ascat)
if ( use_iasi     )  systemParameters % varList_iasi(:)     = varList_iasi(1:varListSize_iasi)
if ( use_oscat    )  systemParameters % varList_oscat(:)    = varList_oscat(1:varListSize_oscat)
if ( use_windsat  )  systemParameters % varList_windsat(:)  = varList_windsat(1:varListSize_windsat)
if ( use_cygnss   )  systemParameters % varList_cygnss(:)   = varList_cygnss(1:varListSize_cygnss)

if ( use_sound    )  systemParameters % use_varList_sound(:)    = use_varList_sound(1:varListSize_sound)
if ( use_airep    )  systemParameters % use_varList_airep(:)    = use_varList_airep(1:varListSize_airep)
if ( use_synop    )  systemParameters % use_varList_synop(:)    = use_varList_synop(1:varListSize_synop)
if ( use_amv      )  systemParameters % use_varList_amv(:)      = use_varList_amv(1:varListSize_amv)
if ( use_gpsro    )  systemParameters % use_varList_gpsro(:)    = use_varList_gpsro(1:varListSize_gpsro)
if ( use_airs     )  systemParameters % use_varList_airs(:)     = use_varList_airs(1:varListSize_airs)
if ( use_quikscat )  systemParameters % use_varList_quikscat(:) = use_varList_quikscat(1:varListSize_quikscat)
if ( use_ascat    )  systemParameters % use_varList_ascat(:)    = use_varList_ascat(1:varListSize_ascat)
if ( use_iasi     )  systemParameters % use_varList_iasi(:)     = use_varList_iasi(1:varListSize_iasi)
if ( use_oscat    )  systemParameters % use_varList_oscat(:)    = use_varList_oscat(1:varListSize_oscat)
if ( use_windsat  )  systemParameters % use_varList_windsat(:)  = use_varList_windsat(1:varListSize_windsat)
if ( use_cygnss   )  systemParameters % use_varList_cygnss(:)   = use_varList_cygnss(1:varListSize_cygnss)
systemParameters % rd = rd
systemParameters % rc = rc
systemParameters % rd_z = rd_z
systemParameters % rc_z = rc_z
systemParameters % inflationFactor = inflationFactor


!================================================
return
stop
end subroutine getSystemParameter
