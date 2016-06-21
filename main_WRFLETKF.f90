
include 'mod_derivedType.f90'
!include 'mod_basicUtility.f90'  ! compiled individually for default real as 8-byte.
include 'mod_systemUtility.f90'
include 'mod_IOUtility.f90'
include 'mod_math.f90'
include 'mod_assimilationUtility.f90'

use derivedType
use basicUtility
use IOUtility
use systemUtility
use assimilationUtility
use omp_lib  ! for OpenMP (internal)

implicit none

type(systemParameter)            :: systemParameters
type(domainInfo),allocatable     :: domain(:)
type(domainInfo)                 :: domain_mean
type(obsParent)                  :: sounding,synop,amv,gpsro
type(backgroundInfo),allocatable :: background(:)
type(backgroundInfo),allocatable :: analysis(:)
integer                          :: ensembleSize

type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachMassGrid => null()
type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachUGrid    => null()
type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachVGrid    => null()
type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachWGrid    => null()

integer io,iwe,isn,iz

real         ct0,ct1
real(kind=8) wt0,wt1,walltime_assimilation
data walltime_assimilation/0d0/
!================================================

print*,'Getting System Parameters...'
call getSystemParameter(systemParameters)
print*,'Done.'
print*,repeat('=',20)

ensembleSize = systemParameters % ensembleSize

allocate( domain(ensembleSize) )

print*,'Getting Domain...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call getDomain(domain(:),ensembleSize)
call build_meanDomain(domain(:),ensembleSize,domain_mean)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(Get domain) =',ct1-ct0,'sec'
print*,'walltime(Get domain) =',wt1-wt0,'sec'


print*,repeat('=',20)
print*,'Getting Observations...'
if ( systemParameters % use_sound )  call getSounding(sounding , systemParameters%varList_sound(:) , systemParameters%varListSize_sound , systemParameters%use_varList_sound )
if ( systemParameters % use_synop )  call getSynop(synop, systemParameters%varList_synop(:) , systemParameters%varListSize_synop , systemParameters%use_varList_synop )
if ( systemParameters % use_amv )    call getAMV(amv, systemParameters%varList_amv(:) , systemParameters%varListSize_amv , systemParameters%use_varList_amv )
if ( systemParameters % use_gpsro )  call getGPSRO(gpsro, systemParameters%varList_gpsro(:) , systemParameters%varListSize_gpsro , systemParameters%use_varList_gpsro )
print*,'Done.'
if ( systemParameters % use_sound )  print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) set to be unavailable by default.'
if ( systemParameters % use_synop )  print*,'There are ',count(.not.synop%obs(:)%available),'/',synop%obsNum,'synop(s) set to be unavailable by default.'
if ( systemParameters % use_amv )    print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) set to be unavailable by default.'
if ( systemParameters % use_gpsro )  print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) set to be unavailable by default.'


wt0 = omp_get_wtime()
call cpu_time(ct0)

print*,repeat('=',20)
print*,'Checking if observations inside horizontal domain...'
if ( systemParameters % use_sound )  call check_ifObsInsideHorizontalDomain(domain(1),sounding)
if ( systemParameters % use_synop )  call check_ifObsInsideHorizontalDomain(domain(1),synop)
if ( systemParameters % use_amv )    call check_ifObsInsideHorizontalDomain(domain(1),amv)
if ( systemParameters % use_gpsro )  call check_ifObsInsideHorizontalDomain(domain(1),gpsro)

if ( systemParameters % use_sound )  print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavailable.'
if ( systemParameters % use_synop )  print*,'There are ',count(.not.synop%obs(:)%available),'/',synop%obsNum,'synop(s) unavailable.'
if ( systemParameters % use_amv )    print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) unavailable.'
if ( systemParameters % use_gpsro )  print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) unavailable.'


print*,repeat('=',20)
print*,'Turning observations with invalid value into unavailable...'
if ( systemParameters % use_sound )  call turnObsWithInvalidValueIntoUnavailable(sounding)
if ( systemParameters % use_amv )    call turnObsWithInvalidValueIntoUnavailable(amv)
if ( systemParameters % use_gpsro )  call turnObsWithInvalidValueIntoUnavailable(gpsro)

if ( systemParameters % use_sound )  print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavailable.'
if ( systemParameters % use_amv )    print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) unavailable.'
if ( systemParameters % use_gpsro )  print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) unavailable.'


print*,repeat('=',20)
print*,'Checking if observations inside vertical domain...'
if ( systemParameters % use_sound )  call check_ifObsInsideVerticalDomain(domain(:),ensembleSize,sounding)
if ( systemParameters % use_amv )    call check_ifObsInsideVerticalDomain(domain(:),ensembleSize,amv)
if ( systemParameters % use_gpsro )  call check_ifObsInsideVerticalDomain(domain(:),ensembleSize,gpsro)

if ( systemParameters % use_sound )  print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavailable.'
if ( systemParameters % use_amv )    print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) unavailable.'
if ( systemParameters % use_gpsro )  print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) unavailable.'

call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'cpu time (Check observations)=',ct1-ct0,'sec'
print*,'walltime (Check observations)=',wt1-wt0,'sec'


print*,repeat('=',20)
wt0 = omp_get_wtime()
allocate( background(ensembleSize) )
call getBackground(background(:),ensembleSize,domain(:))
wt1 = omp_get_wtime()
print*,'walltime(get background) =',wt1-wt0,'sec'

print*,repeat('=',20)
print*,'Converting background to tempertature...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call convertBackgroundToTemperature(background(:),ensembleSize,domain(:))
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(background to T) =',ct1-ct0,'sec'
print*,'walltime(background to T) =',wt1-wt0,'sec'


print*,repeat('=',20)
print*,'Converting background to sounding...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call convertBackgroundToSounding(background(:),ensembleSize,domain(:),sounding)
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


print*,repeat('=',20)
print*,'Rearranging observation base on their veritcal position...'
wt0 = omp_get_wtime()
call rearrangeObsToBeBottomToTop(sounding)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'walltime(rearrange observation) =',wt1-wt0,'sec'


print*,repeat('=',20)
print*,'Initilizing the analysis...'
wt0 = omp_get_wtime()
allocate( analysis(ensembleSize) )
call initializeAnalysis(background,analysis,ensembleSize,domain)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'walltime(initialize analysis) =',wt1-wt0,'sec'


!
!  Assimilation on mass grid.
!
print*,repeat('=',20)
print*,'Mapping observations to each mass grid...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call mapObsToEachMassGrid(obsListOfEachMassGrid,sounding,domain_mean,systemParameters)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time (mapObsToEachMassGrid)=',ct1-ct0,'sec'
print*,'walltime (mapObsToEachMassGrid)=',wt1-wt0,'sec'
print*,'There are',count(obsListOfEachMassGrid(:,:,:)%vectorSize>0),'grids may need to be updated.'
print*,'Max obs per grid=',maxval(obsListOfEachMassGrid(:,:,:)%vectorSize)
print*,'Total obs for all grids=',sum(obsListOfEachMassGrid(:,:,:)%vectorSize)


print*,'Starting assimilation...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call assimilate_massGrid(background(:),analysis(:),ensembleSize,domain(:),domain_mean,sounding,obsListOfEachMassGrid,systemParameters)
call cpu_time(ct1)
wt1 = omp_get_wtime()
walltime_assimilation = walltime_assimilation + (wt1-wt0)
print*,'Done.'
print*,'cpu time(assimilation on mass grid) =',ct1-ct0,'sec'
print*,'walltime(assimilation on mass grid) =',wt1-wt0,'sec'

call deallocate_obsListOfEachGrid(obsListOfEachMassGrid)


!
!  Assimilation on u grid.
!
print*,repeat('=',20)
print*,'Mapping observations to each u grid...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call mapObsToEachUGrid(obsListOfEachUGrid,sounding,domain_mean,systemParameters)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time (mapObsToEachUGrid)=',ct1-ct0,'sec'
print*,'walltime (mapObsToEachUGrid)=',wt1-wt0,'sec'
print*,'There are',count(obsListOfEachUGrid(:,:,:)%vectorSize>0),'grids may need to be updated.'
print*,'Max obs per grid=',maxval(obsListOfEachUGrid(:,:,:)%vectorSize)
print*,'Total obs for all grids=',sum(obsListOfEachUGrid(:,:,:)%vectorSize)


print*,'Starting assimilation...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call assimilate_uGrid(background(:),analysis(:),ensembleSize,domain(:),domain_mean,sounding,obsListOfEachUGrid,systemParameters)
call cpu_time(ct1)
wt1 = omp_get_wtime()
walltime_assimilation = walltime_assimilation + (wt1-wt0)
print*,'Done.'
print*,'cpu time(assimilation on u grid) =',ct1-ct0,'sec'
print*,'walltime(assimilation on u grid) =',wt1-wt0,'sec'

call deallocate_obsListOfEachGrid(obsListOfEachUGrid)


!
!  Assimilation on v grid.
!
print*,repeat('=',20)
print*,'Mapping observations to each v grid...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call mapObsToEachVGrid(obsListOfEachVGrid,sounding,domain_mean,systemParameters)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time (mapObsToEachVGrid)=',ct1-ct0,'sec'
print*,'walltime (mapObsToEachVGrid)=',wt1-wt0,'sec'
print*,'There are',count(obsListOfEachVGrid(:,:,:)%vectorSize>0),'grids may need to be updated.'
print*,'Max obs per grid=',maxval(obsListOfEachVGrid(:,:,:)%vectorSize)
print*,'Total obs for all grids=',sum(obsListOfEachVGrid(:,:,:)%vectorSize)


print*,'Starting assimilation...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call assimilate_vGrid(background(:),analysis(:),ensembleSize,domain(:),domain_mean,sounding,obsListOfEachVGrid,systemParameters)
call cpu_time(ct1)
wt1 = omp_get_wtime()
walltime_assimilation = walltime_assimilation + (wt1-wt0)
print*,'Done.'
print*,'cpu time(assimilation on v grid) =',ct1-ct0,'sec'
print*,'walltime(assimilation on v grid) =',wt1-wt0,'sec'

call deallocate_obsListOfEachGrid(obsListOfEachVGrid)


!
!  Assimilation on w grid.
!
print*,repeat('=',20)
print*,'Mapping observations to each w grid...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call mapObsToEachWGrid(obsListOfEachWGrid,sounding,domain_mean,systemParameters)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time (mapObsToEachWGrid)=',ct1-ct0,'sec'
print*,'walltime (mapObsToEachWGrid)=',wt1-wt0,'sec'
print*,'There are',count(obsListOfEachWGrid(:,:,:)%vectorSize>0),'grids may need to be updated.'
print*,'Max obs per grid=',maxval(obsListOfEachWGrid(:,:,:)%vectorSize)
print*,'Total obs for all grids=',sum(obsListOfEachWGrid(:,:,:)%vectorSize)


print*,'Starting assimilation...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call assimilate_wGrid(background(:),analysis(:),ensembleSize,domain(:),domain_mean,sounding,obsListOfEachWGrid,systemParameters)
call cpu_time(ct1)
wt1 = omp_get_wtime()
walltime_assimilation = walltime_assimilation + (wt1-wt0)
print*,'Done.'
print*,'cpu time(assimilation on w grid) =',ct1-ct0,'sec'
print*,'walltime(assimilation on w grid) =',wt1-wt0,'sec'

call deallocate_obsListOfEachGrid(obsListOfEachWGrid)


print*,repeat('=',20)
print*,'walltime(all assimilation process) =',walltime_assimilation,'sec'
print*,repeat('=',20)
!
!  Output analysis.
!
wt0 = omp_get_wtime()
call cpu_time(ct0)
call outputAnalysis(analysis(:),ensembleSize,domain(:))
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'cpu time(outputAnalysis) =',ct1-ct0,'sec'
print*,'walltime(outputAnalysis) =',wt1-wt0,'sec'



!================================================
stop
end
