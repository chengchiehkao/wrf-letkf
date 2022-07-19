
use derivedType
use IOUtility
use systemUtility
use assimilationUtility
use omp_lib  ! OpenMP's library

implicit none

type(systemParameter)            :: systemParameters
type(domainInfo),allocatable     :: domain(:,:)
type(domainInfo),allocatable     :: domain_mean(:)
type(obsParent)                  :: sounding,airep,synop,amv,gpsro,airs,quikscat,ascat,iasi,oscat,windsat,cygnss
type(obsParent)                  :: allObs
type(backgroundInfo),allocatable :: background(:,:)
type(backgroundInfo),allocatable :: analysis(:,:)
integer                          :: ensembleSize

type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachMassGrid => null()
type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachUGrid    => null()
type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachVGrid    => null()
type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachWGrid    => null()

integer io,iwe,isn,iz,domainID

real         ct0,ct1
real(kind=8) wt0,wt1,walltime_assimilation
data walltime_assimilation/0d0/
!================================================

print*,'Getting System Parameters...'
call getSystemParameter(systemParameters)
print*,'Done.'
print*,repeat('=',20)

ensembleSize = systemParameters % ensembleSize

allocate( domain(ensembleSize,systemParameters%max_domain) )
allocate( domain_mean(systemParameters%max_domain) )

print*,'Getting Domain...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
do domainID = 1 , systemParameters%max_domain
    call getDomain(domain(:,domainID),ensembleSize,domainID)
    call build_pressureU(domain(:,domainID),ensembleSize)
    call build_pressureV(domain(:,domainID),ensembleSize)
    call convertGPHToGMH(domain(:,domainID),ensembleSize)
    call build_meanDomain(domain(:,domainID),ensembleSize,domain_mean(domainID))
enddo
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(Get domain) =',ct1-ct0,'sec'
print*,'walltime(Get domain) =',wt1-wt0,'sec'


print*,repeat('=',20)
print*,'Getting Observations...'
if ( systemParameters % use_sound    )  call getSounding(sounding , systemParameters)
if ( systemParameters % use_airep    )  call getAirep(airep , systemParameters)
if ( systemParameters % use_synop    )  call getSynop(synop , systemParameters)
if ( systemParameters % use_amv      )  call getAMV(amv , systemParameters)
if ( systemParameters % use_gpsro    )  call getGPSRO(gpsro , systemParameters)
if ( systemParameters % use_airs     )  call getAIRS(airs , systemParameters)
if ( systemParameters % use_quikscat )  call getQuikSCAT(quikscat , systemParameters)
if ( systemParameters % use_ascat    )  call getASCAT(ascat , systemParameters)
if ( systemParameters % use_iasi     )  call getIASI(iasi , systemParameters)
if ( systemParameters % use_oscat    )  call getOSCAT(oscat , systemParameters)
if ( systemParameters % use_windsat  )  call getWindSat(windsat , systemParameters)
if ( systemParameters % use_cygnss   )  call getCYGNSS(cygnss , systemParameters)
print*,'Done.'
if ( systemParameters % use_sound    )  print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) set to be unavailable by default.'
if ( systemParameters % use_airep    )  print*,'There are ',count(.not.airep%obs(:)%available),'/',airep%obsNum,'airep(s) set to be unavailable by default.'
if ( systemParameters % use_synop    )  print*,'There are ',count(.not.synop%obs(:)%available),'/',synop%obsNum,'synop(s) set to be unavailable by default.'
if ( systemParameters % use_amv      )  print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) set to be unavailable by default.'
if ( systemParameters % use_gpsro    )  print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) set to be unavailable by default.'
if ( systemParameters % use_airs     )  print*,'There are ',count(.not.airs%obs(:)%available),'/',airs%obsNum,'airs(s) set to be unavailable by default.'
if ( systemParameters % use_quikscat )  print*,'There are ',count(.not.quikscat%obs(:)%available),'/',quikscat%obsNum,'quikscat(s) set to be unavailable by default.'
if ( systemParameters % use_ascat    )  print*,'There are ',count(.not.ascat%obs(:)%available),'/',ascat%obsNum,'ascat(s) set to be unavailable by default.'
if ( systemParameters % use_iasi     )  print*,'There are ',count(.not.iasi%obs(:)%available),'/',iasi%obsNum,'iasi(s) set to be unavailable by default.'
if ( systemParameters % use_oscat    )  print*,'There are ',count(.not.oscat%obs(:)%available),'/',oscat%obsNum,'oscat(s) set to be unavailable by default.'
if ( systemParameters % use_windsat  )  print*,'There are ',count(.not.windsat%obs(:)%available),'/',windsat%obsNum,'windsat(s) set to be unavailable by default.'
if ( systemParameters % use_cygnss   )  print*,'There are ',count(.not.cygnss%obs(:)%available),'/',cygnss%obsNum,'cygnss(s) set to be unavailable by default.'


wt0 = omp_get_wtime()
call cpu_time(ct0)

print*,repeat('=',20)
print*,'Checking if observations inside horizontal domain...'
do domainID = 1 , systemParameters%max_domain
    if ( systemParameters % use_sound    )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),sounding,domainID)
    if ( systemParameters % use_airep    )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),airep   ,domainID)
    if ( systemParameters % use_synop    )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),synop   ,domainID)
    if ( systemParameters % use_amv      )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),amv     ,domainID)
    if ( systemParameters % use_gpsro    )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),gpsro   ,domainID)
    if ( systemParameters % use_airs     )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),airs    ,domainID)
    if ( systemParameters % use_quikscat )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),quikscat,domainID)
    if ( systemParameters % use_ascat    )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),ascat   ,domainID)
    if ( systemParameters % use_iasi     )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),iasi    ,domainID)
    if ( systemParameters % use_oscat    )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),oscat   ,domainID)
    if ( systemParameters % use_windsat  )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),windsat ,domainID)
    if ( systemParameters % use_cygnss   )  call check_ifObsInsideHorizontalDomain(domain(1,domainID),cygnss  ,domainID)
enddo

if ( systemParameters % use_sound    )  print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavailable.'
if ( systemParameters % use_airep    )  print*,'There are ',count(.not.airep%obs(:)%available),'/',airep%obsNum,'airep(s) unavailable.'
if ( systemParameters % use_synop    )  print*,'There are ',count(.not.synop%obs(:)%available),'/',synop%obsNum,'synop(s) unavailable.'
if ( systemParameters % use_amv      )  print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) unavailable.'
if ( systemParameters % use_gpsro    )  print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) unavailable.'
if ( systemParameters % use_airs     )  print*,'There are ',count(.not.airs%obs(:)%available),'/',airs%obsNum,'airs(s) unavailable.'
if ( systemParameters % use_quikscat )  print*,'There are ',count(.not.quikscat%obs(:)%available),'/',quikscat%obsNum,'quikscat(s) unavailable.'
if ( systemParameters % use_ascat    )  print*,'There are ',count(.not.ascat%obs(:)%available),'/',ascat%obsNum,'ascat(s) unavailable.'
if ( systemParameters % use_iasi     )  print*,'There are ',count(.not.iasi%obs(:)%available),'/',iasi%obsNum,'iasi(s) unavailable.'
if ( systemParameters % use_oscat    )  print*,'There are ',count(.not.oscat%obs(:)%available),'/',oscat%obsNum,'oscat(s) unavailable.'
if ( systemParameters % use_windsat  )  print*,'There are ',count(.not.windsat%obs(:)%available),'/',windsat%obsNum,'windsat(s) unavailable.'
if ( systemParameters % use_cygnss   )  print*,'There are ',count(.not.cygnss%obs(:)%available),'/',cygnss%obsNum,'cygnss(s) unavailable.'


print*,repeat('=',20)
print*,'Turning observations with invalid value into unavailable...'
if ( systemParameters % use_sound    )  call turnObsWithInvalidValueIntoUnavailable(sounding)
if ( systemParameters % use_airep    )  call turnObsWithInvalidValueIntoUnavailable(airep)
if ( systemParameters % use_amv      )  call turnObsWithInvalidValueIntoUnavailable(amv)
if ( systemParameters % use_gpsro    )  call turnObsWithInvalidValueIntoUnavailable(gpsro)
if ( systemParameters % use_airs     )  call turnObsWithInvalidValueIntoUnavailable(airs)
if ( systemParameters % use_quikscat )  call turnObsWithInvalidValueIntoUnavailable(quikscat)
if ( systemParameters % use_ascat    )  call turnObsWithInvalidValueIntoUnavailable(ascat)
if ( systemParameters % use_iasi     )  call turnObsWithInvalidValueIntoUnavailable(iasi)
if ( systemParameters % use_oscat    )  call turnObsWithInvalidValueIntoUnavailable(oscat)
if ( systemParameters % use_windsat  )  call turnObsWithInvalidValueIntoUnavailable(windsat)
if ( systemParameters % use_cygnss   )  call turnObsWithInvalidValueIntoUnavailable(cygnss)

if ( systemParameters % use_sound    )  print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavailable.'
if ( systemParameters % use_airep    )  print*,'There are ',count(.not.airep%obs(:)%available),'/',airep%obsNum,'airep(s) unavailable.'
if ( systemParameters % use_amv      )  print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) unavailable.'
if ( systemParameters % use_gpsro    )  print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) unavailable.'
if ( systemParameters % use_airs     )  print*,'There are ',count(.not.airs%obs(:)%available),'/',airs%obsNum,'airs(s) unavailable.'
if ( systemParameters % use_quikscat )  print*,'There are ',count(.not.quikscat%obs(:)%available),'/',quikscat%obsNum,'quikscat(s) unavailable.'
if ( systemParameters % use_ascat    )  print*,'There are ',count(.not.ascat%obs(:)%available),'/',ascat%obsNum,'ascat(s) unavailable.'
if ( systemParameters % use_iasi     )  print*,'There are ',count(.not.iasi%obs(:)%available),'/',iasi%obsNum,'iasi(s) unavailable.'
if ( systemParameters % use_oscat    )  print*,'There are ',count(.not.oscat%obs(:)%available),'/',oscat%obsNum,'oscat(s) unavailable.'
if ( systemParameters % use_windsat  )  print*,'There are ',count(.not.windsat%obs(:)%available),'/',windsat%obsNum,'windsat(s) unavailable.'
if ( systemParameters % use_cygnss   )  print*,'There are ',count(.not.cygnss%obs(:)%available),'/',cygnss%obsNum,'cygnss(s) unavailable.'


print*,repeat('=',20)
print*,'Checking if observations inside vertical domain...'
do domainID = 1 , systemParameters%max_domain
    if ( systemParameters % use_sound )  call check_ifObsInsideVerticalDomain(domain(:,domainID),ensembleSize,sounding,domainID)
    if ( systemParameters % use_airep )  call check_ifObsInsideVerticalDomain(domain(:,domainID),ensembleSize,airep   ,domainID)
    if ( systemParameters % use_amv )    call check_ifObsInsideVerticalDomain(domain(:,domainID),ensembleSize,amv     ,domainID)
    if ( systemParameters % use_gpsro )  call check_ifObsInsideVerticalDomain(domain(:,domainID),ensembleSize,gpsro   ,domainID)
    if ( systemParameters % use_airs )   call check_ifObsInsideVerticalDomain(domain(:,domainID),ensembleSize,airs    ,domainID)
    if ( systemParameters % use_iasi )   call check_ifObsInsideVerticalDomain(domain(:,domainID),ensembleSize,iasi    ,domainID)
    !  QuikSCAT does NOT have to do vertical check because it's always 10-m higher than surface.
    !  ASCAT    does NOT have to do vertical check because it's always 10-m higher than surface.
    !  OSCAT    does NOT have to do vertical check because it's always 10-m higher than surface.
    !  WindSat  does NOT have to do vertical check because it's always 10-m higher than surface.
    !  CYGNSS   does NOT have to do vertical check because it's always 10-m higher than surfcae.
enddo

if ( systemParameters % use_sound    )  print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavailable.'
if ( systemParameters % use_airep    )  print*,'There are ',count(.not.airep%obs(:)%available),'/',airep%obsNum,'airep(s) unavailable.'
if ( systemParameters % use_amv      )  print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) unavailable.'
if ( systemParameters % use_gpsro    )  print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) unavailable.'
if ( systemParameters % use_airs     )  print*,'There are ',count(.not.airs%obs(:)%available),'/',airs%obsNum,'airs(s) unavailable.'
if ( systemParameters % use_quikscat )  print*,'There are ',count(.not.quikscat%obs(:)%available),'/',quikscat%obsNum,'quikscat(s) unavailable.'
if ( systemParameters % use_ascat    )  print*,'There are ',count(.not.ascat%obs(:)%available),'/',ascat%obsNum,'ascat(s) unavailable.'
if ( systemParameters % use_iasi     )  print*,'There are ',count(.not.iasi%obs(:)%available),'/',iasi%obsNum,'iasi(s) unavailable.'
if ( systemParameters % use_oscat    )  print*,'There are ',count(.not.oscat%obs(:)%available),'/',oscat%obsNum,'oscat(s) unavailable.'
if ( systemParameters % use_windsat  )  print*,'There are ',count(.not.windsat%obs(:)%available),'/',windsat%obsNum,'windsat(s) unavailable.'
if ( systemParameters % use_cygnss    )  print*,'There are ',count(.not.cygnss%obs(:)%available),'/',cygnss%obsNum,'cygnss(s) unavailable.'

call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'cpu time (Check observations)=',ct1-ct0,'sec'
print*,'walltime (Check observations)=',wt1-wt0,'sec'


print*,repeat('=',20)
wt0 = omp_get_wtime()
allocate( background(ensembleSize,systemParameters%max_domain) )
do domainID = 1 , systemParameters%max_domain
    call getBackground(background(:,domainID),ensembleSize,domain(:,domainID),domainID)
enddo
wt1 = omp_get_wtime()
print*,'walltime(get background) =',wt1-wt0,'sec'

print*,repeat('=',20)
print*,'Converting background to tempertature...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
do domainID = 1 , systemParameters%max_domain
    call convertBackgroundToTemperature(background(:,domainID),ensembleSize,domain(:,domainID))
enddo
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(background to T) =',ct1-ct0,'sec'
print*,'walltime(background to T) =',wt1-wt0,'sec'

print*,repeat('=',20)
print*,'Converting background to relative humidity...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
do domainID = 1 , systemParameters%max_domain
     call convertBackgroundToRH(background(:,domainID),ensembleSize,domain(:,domainID))
enddo
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(background to RH) =',ct1-ct0,'sec'
print*,'walltime(background to RH) =',wt1-wt0,'sec'


if ( systemParameters % use_sound ) then
    print*,repeat('=',20)
    print*,'Converting background to sounding...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToSounding(background(:,domainID),ensembleSize,domain(:,domainID),sounding,systemParameters,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of sounding) =',ct1-ct0,'sec'
    print*,'walltime(H of sounding) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of sounding...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setSoundingError(sounding)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set sounding error) =',ct1-ct0,'sec'
    print*,'walltime(set sounding error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavailable.'
endif


if ( systemParameters % use_airep ) then
    print*,repeat('=',20)
    print*,'Converting background to airep...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToAirep(background(:,domainID),ensembleSize,domain(:,domainID),airep,systemParameters,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of airep) =',ct1-ct0,'sec'
    print*,'walltime(H of airep) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of airep...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setAirepError(airep)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set airep error) =',ct1-ct0,'sec'
    print*,'walltime(set airep error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.airep%obs(:)%available),'/',airep%obsNum,'airep(s) unavailable.'
endif


if ( systemParameters % use_amv ) then
    print*,repeat('=',20)
    print*,'Converting background to AMV...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToAMV(background(:,domainID),ensembleSize,domain(:,domainID),amv,systemParameters,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of AMV) =',ct1-ct0,'sec'
    print*,'walltime(H of AMV) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of AMV...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setAMVError(amv)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set AMV error) =',ct1-ct0,'sec'
    print*,'walltime(set AMV error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'AMV(s) unavailable.'
endif


if ( systemParameters % use_gpsro ) then
    print*,repeat('=',20)
    print*,'Converting background to GPSRO...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToRefractivity(background(:,domainID),ensembleSize,domain(:,domainID))
        call convertBackgroundToGPSRO(background(:,domainID),ensembleSize,domain(:,domainID),gpsro,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of GPSRO) =',ct1-ct0,'sec'
    print*,'walltime(H of GPSRO) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of GPSRO...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setGPSROError(gpsro)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set GPSRO error) =',ct1-ct0,'sec'
    print*,'walltime(set GPSRO error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'GPSRO(s) unavailable.'

    print*,repeat('=',20)
    print*,'Mapping GPSRO from GMH to P...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call mapGPSROFromGMHToP(domain_mean(domainID),gpsro,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(map GPSRO) =',ct1-ct0,'sec'
    print*,'walltime(map GPSRO) =',wt1-wt0,'sec'
endif



if ( systemParameters % use_airs ) then
    print*,repeat('=',20)
    print*,'Converting background to AIRS...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToAIRS(background(:,domainID),ensembleSize,domain(:,domainID),airs,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of AIRS) =',ct1-ct0,'sec'
    print*,'walltime(H of AIRS) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of AIRS...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setAIRSError(airs)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set AIRS error) =',ct1-ct0,'sec'
    print*,'walltime(set AIRS error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.airs%obs(:)%available),'/',airs%obsNum,'AIRS(s) unavailable.'
endif


if ( systemParameters % use_quikscat ) then
    print*,repeat('=',20)
    print*,'Converting background to QuikSCAT...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToQuikSCAT(background(:,domainID),ensembleSize,domain(:,domainID),domain_mean(domainID),quikscat,systemParameters,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of QuikSCAT) =',ct1-ct0,'sec'
    print*,'walltime(H of QuikSCAT) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of QuikSCAT...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setQuikSCATError(quikscat)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set QuikSCAT error) =',ct1-ct0,'sec'
    print*,'walltime(set QuikSCAT error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.quikscat%obs(:)%available),'/',quikscat%obsNum,'QuikSCAT(s) unavailable.'
endif


if ( systemParameters % use_ascat ) then
    print*,repeat('=',20)
    print*,'Converting background to ASCAT...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToASCAT(background(:,domainID),ensembleSize,domain(:,domainID),domain_mean(domainID),ascat,systemParameters,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of ASCAT) =',ct1-ct0,'sec'
    print*,'walltime(H of ASCAT) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of ASCAT...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setASCATError(ascat)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set ASCAT error) =',ct1-ct0,'sec'
    print*,'walltime(set ASCAT error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.ascat%obs(:)%available),'/',ascat%obsNum,'ASCAT(s) unavailable.'
endif


if ( systemParameters % use_iasi ) then
    print*,repeat('=',20)
    print*,'Converting background to IASI...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToIASI(background(:,domainID),ensembleSize,domain(:,domainID),iasi,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of IASI) =',ct1-ct0,'sec'
    print*,'walltime(H of IASI) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of IASI...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setIASIError(iasi)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set IASI error) =',ct1-ct0,'sec'
    print*,'walltime(set IASI error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.iasi%obs(:)%available),'/',iasi%obsNum,'IASI(s) unavailable.'
endif


if ( systemParameters % use_oscat ) then
    print*,repeat('=',20)
    print*,'Converting background to OSCAT...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToOSCAT(background(:,domainID),ensembleSize,domain(:,domainID),domain_mean(domainID),oscat,systemParameters,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of OSCAT) =',ct1-ct0,'sec'
    print*,'walltime(H of OSCAT) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of OSCAT...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setOSCATError(oscat)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set OSCAT error) =',ct1-ct0,'sec'
    print*,'walltime(set OSCAT error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.oscat%obs(:)%available),'/',oscat%obsNum,'OSCAT(s) unavailable.'
endif


if ( systemParameters % use_windsat ) then
    print*,repeat('=',20)
    print*,'Converting background to WindSat...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToWindSat(background(:,domainID),ensembleSize,domain(:,domainID),domain_mean(domainID),windsat,systemParameters,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of WindSat) =',ct1-ct0,'sec'
    print*,'walltime(H of WindSat) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of WindSat...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setWindSatError(windsat)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set WindSat error) =',ct1-ct0,'sec'
    print*,'walltime(set WindSat error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.windsat%obs(:)%available),'/',windsat%obsNum,'WindSat(s) unavailable.'
endif


if ( systemParameters % use_cygnss ) then
    print*,repeat('=',20)
    print*,'Converting background to CYGNSS...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    do domainID = 1,systemParameters%max_domain
        call convertBackgroundToCYGNSS(background(:,domainID),ensembleSize,domain(:,domainID),domain_mean(domainID),cygnss,domainID)
    enddo
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(H of CYGNSS) =',ct1-ct0,'sec'
    print*,'walltime(H of CYGNSS) =',wt1-wt0,'sec'


    print*,repeat('=',20)
    print*,'Setting error of CYGNSS...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call setCYGNSSError(cygnss)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time(set CYGNSS error) =',ct1-ct0,'sec'
    print*,'walltime(set CYGNSS error) =',wt1-wt0,'sec'
    print*,'There are ',count(.not.cygnss%obs(:)%available),'/',cygnss%obsNum,'CYGNSS(s) unavailable.'
endif


print*,repeat('=',20)
print*,'Merging observation(s)...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
if ( systemParameters % use_sound    )  call mergeObs(allObs,sounding)
if ( systemParameters % use_airep    )  call mergeObs(allObs,airep)
if ( systemParameters % use_amv      )  call mergeObs(allObs,amv)
if ( systemParameters % use_gpsro    )  call mergeObs(allObs,gpsro)
if ( systemParameters % use_airs     )  call mergeObs(allObs,airs)
if ( systemParameters % use_quikscat )  call mergeObs(allObs,quikscat)
if ( systemParameters % use_ascat    )  call mergeObs(allObs,ascat)
if ( systemParameters % use_iasi     )  call mergeObs(allObs,iasi)
if ( systemParameters % use_oscat    )  call mergeObs(allObs,oscat)
if ( systemParameters % use_windsat  )  call mergeObs(allObs,windsat)
if ( systemParameters % use_cygnss   )  call mergeObs(allObs,cygnss)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(merge obs) =',ct1-ct0,'sec'
print*,'walltime(merge obs) =',wt1-wt0,'sec'
print*,'There are ',count(.not.allObs%obs(:)%available),'/',allObs%obsNum,'obs unavailable.'


print*,repeat('=',20)
print*,'Rearranging observation base on their veritcal position...'
wt0 = omp_get_wtime()
call rearrangeObsToBeBottomToTop(allObs)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'walltime(rearrange observation) =',wt1-wt0,'sec'


print*,repeat('=',20)
print*,'Initilizing the analysis...'
wt0 = omp_get_wtime()
allocate( analysis(ensembleSize,systemParameters%max_domain) )
do domainID = 1,systemParameters%max_domain 
    call initializeAnalysis(background(:,domainID),analysis(:,domainID),ensembleSize,domain(:,domainID))
enddo
wt1 = omp_get_wtime()
print*,'Done.'
print*,'walltime(initialize analysis) =',wt1-wt0,'sec'


!
!  Assimilation on mass grid.
!

print*,repeat('=',20)
do domainID = 1,systemParameters%max_domain
    write(*,'("D",i2.2,": ")',advance='no' ) domainID
    print*,'Mapping observations to each mass grid...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call mapObsToEachMassGrid(obsListOfEachMassGrid,allObs,domain_mean(domainID),systemParameters)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time (mapObsToEachMassGrid)=',ct1-ct0,'sec'
    print*,'walltime (mapObsToEachMassGrid)=',wt1-wt0,'sec'
    print*,'There are',count(obsListOfEachMassGrid(:,:,:)%vectorSize>0),'grids may need to be updated.'
    print*,'Max obs per grid=',maxval(obsListOfEachMassGrid(:,:,:)%vectorSize)
    print*,'Total obs for all grids=',sum(obsListOfEachMassGrid(:,:,:)%vectorSize)


    write(*,'("D",i2.2,": ")',advance='no' ) domainID
    print*,'Starting assimilation...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call assimilate_massGrid(background(:,domainID),analysis(:,domainID),ensembleSize,domain_mean(domainID),allObs,obsListOfEachMassGrid,systemParameters,domainID)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    walltime_assimilation = walltime_assimilation + (wt1-wt0)
    print*,'Done.'
    print*,'cpu time(assimilation on mass grid) =',ct1-ct0,'sec'
    print*,'walltime(assimilation on mass grid) =',wt1-wt0,'sec'

    call deallocate_obsListOfEachGrid(obsListOfEachMassGrid)
enddo


!
!  Assimilation on u grid.
!
print*,repeat('=',20)
do domainID = 1,systemParameters%max_domain
    write(*,'("D",i2.2,": ")',advance='no' ) domainID
    print*,'Mapping observations to each u grid...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call mapObsToEachUGrid(obsListOfEachUGrid,allObs,domain_mean(domainID),systemParameters)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time (mapObsToEachUGrid)=',ct1-ct0,'sec'
    print*,'walltime (mapObsToEachUGrid)=',wt1-wt0,'sec'
    print*,'There are',count(obsListOfEachUGrid(:,:,:)%vectorSize>0),'grids may need to be updated.'
    print*,'Max obs per grid=',maxval(obsListOfEachUGrid(:,:,:)%vectorSize)
    print*,'Total obs for all grids=',sum(obsListOfEachUGrid(:,:,:)%vectorSize)


    write(*,'("D",i2.2,": ")',advance='no' ) domainID
    print*,'Starting assimilation...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call assimilate_uGrid(background(:,domainID),analysis(:,domainID),ensembleSize,domain_mean(domainID),allObs,obsListOfEachUGrid,systemParameters,domainID)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    walltime_assimilation = walltime_assimilation + (wt1-wt0)
    print*,'Done.'
    print*,'cpu time(assimilation on u grid) =',ct1-ct0,'sec'
    print*,'walltime(assimilation on u grid) =',wt1-wt0,'sec'

    call deallocate_obsListOfEachGrid(obsListOfEachUGrid)
enddo


!
!  Assimilation on v grid.
!
print*,repeat('=',20)
do domainID = 1,systemParameters%max_domain
    write(*,'("D",i2.2,": ")',advance='no' ) domainID
    print*,'Mapping observations to each v grid...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call mapObsToEachVGrid(obsListOfEachVGrid,allObs,domain_mean(domainID),systemParameters)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time (mapObsToEachVGrid)=',ct1-ct0,'sec'
    print*,'walltime (mapObsToEachVGrid)=',wt1-wt0,'sec'
    print*,'There are',count(obsListOfEachVGrid(:,:,:)%vectorSize>0),'grids may need to be updated.'
    print*,'Max obs per grid=',maxval(obsListOfEachVGrid(:,:,:)%vectorSize)
    print*,'Total obs for all grids=',sum(obsListOfEachVGrid(:,:,:)%vectorSize)


    write(*,'("D",i2.2,": ")',advance='no' ) domainID
    print*,'Starting assimilation...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call assimilate_vGrid(background(:,domainID),analysis(:,domainID),ensembleSize,domain_mean(domainID),allObs,obsListOfEachVGrid,systemParameters,domainID)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    walltime_assimilation = walltime_assimilation + (wt1-wt0)
    print*,'Done.'
    print*,'cpu time(assimilation on v grid) =',ct1-ct0,'sec'
    print*,'walltime(assimilation on v grid) =',wt1-wt0,'sec'

    call deallocate_obsListOfEachGrid(obsListOfEachVGrid)
enddo


!
!  Assimilation on w grid.
!
print*,repeat('=',20)
do domainID = 1,systemParameters%max_domain
    write(*,'("D",i2.2,": ")',advance='no' ) domainID
    print*,'Mapping observations to each w grid...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call mapObsToEachWGrid(obsListOfEachWGrid,allObs,domain_mean(domainID),systemParameters)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    print*,'Done.'
    print*,'cpu time (mapObsToEachWGrid)=',ct1-ct0,'sec'
    print*,'walltime (mapObsToEachWGrid)=',wt1-wt0,'sec'
    print*,'There are',count(obsListOfEachWGrid(:,:,:)%vectorSize>0),'grids may need to be updated.'
    print*,'Max obs per grid=',maxval(obsListOfEachWGrid(:,:,:)%vectorSize)
    print*,'Total obs for all grids=',sum(obsListOfEachWGrid(:,:,:)%vectorSize)


    write(*,'("D",i2.2,": ")',advance='no' ) domainID
    print*,'Starting assimilation...'
    wt0 = omp_get_wtime()
    call cpu_time(ct0)
    call assimilate_wGrid(background(:,domainID),analysis(:,domainID),ensembleSize,domain_mean(domainID),allObs,obsListOfEachWGrid,systemParameters,domainID)
    call cpu_time(ct1)
    wt1 = omp_get_wtime()
    walltime_assimilation = walltime_assimilation + (wt1-wt0)
    print*,'Done.'
    print*,'cpu time(assimilation on w grid) =',ct1-ct0,'sec'
    print*,'walltime(assimilation on w grid) =',wt1-wt0,'sec'
    
    call deallocate_obsListOfEachGrid(obsListOfEachWGrid)
enddo


print*,repeat('=',20)
print*,'walltime(all assimilation process) =',walltime_assimilation,'sec'
print*,repeat('=',20)
!
!  Output analysis.
!
wt0 = omp_get_wtime()
call cpu_time(ct0)
do domainID = 1,systemParameters%max_domain
    call outputAnalysis(analysis(:,domainID),ensembleSize,domain(:,domainID),domainID)
enddo
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'cpu time(outputAnalysis) =',ct1-ct0,'sec'
print*,'walltime(outputAnalysis) =',wt1-wt0,'sec'



!================================================
stop
end
