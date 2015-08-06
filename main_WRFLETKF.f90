
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
real(kind=8) wt0,wt1
!================================================

ensembleSize=36

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


print*,'Getting Observations...'
call getSounding(sounding , (/'U         ','V         ','T         ','QVAPOR    '/) , 4 )
call getSynop(synop,        (/'PSFC      ','U         ','V         ','T         ','Td        '/) , 5 )
call getAMV(amv,            (/'U         ','V         '/) , 2 )
call getGPSRO(gpsro,        (/'REF       ','BANGLE    ','ImpctParam'/) , 3 )
sounding%obs(:)%z = 100.d0 * sounding%obs(:)%z  ! convert hPa to Pa
amv%obs(:)%z      = 100.d0 * amv%obs(:)%z       ! convert hPa to Pa
print*,'Done.'

call cpu_time(ct0)

print*,'Checking if observations inside horizontal domain...'
call check_ifObsInsideHorizontalDomain(domain(1),sounding)
call check_ifObsInsideHorizontalDomain(domain(1),synop)
call check_ifObsInsideHorizontalDomain(domain(1),amv)
call check_ifObsInsideHorizontalDomain(domain(1),gpsro)

print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) out of horizontal domain.'
print*,'There are ',count(.not.synop%obs(:)%available),'/',synop%obsNum,'synop(s) out of horizontal domain.'
print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) out of horizontal domain.'
print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) out of horizontal domain.'


print*,repeat('=',20)
print*,'Turning observations with invalid value into unavailable...'
call turnObsWithInvalidValueIntoUnavailable(sounding)
call turnObsWithInvalidValueIntoUnavailable(amv)
call turnObsWithInvalidValueIntoUnavailable(gpsro)

print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavailable.'
print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) unavailable.'
print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) unavailable.'


print*,repeat('=',20)
print*,'Checking if observations inside vertical domain...'
call check_ifObsInsideVerticalDomain(domain(:),ensembleSize,sounding)
call check_ifObsInsideVerticalDomain(domain(:),ensembleSize,amv)
call check_ifObsInsideVerticalDomain(domain(:),ensembleSize,gpsro)

print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavailable.'
print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) unavailable.'
print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) unavailable.'


call cpu_time(ct1)
print*,'cpu time (Check observations)=',ct1-ct0,'sec'



allocate( background(ensembleSize) )
call getBackground(background(:),ensembleSize,domain(:))

print*,'Converting background to tempertature...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call convertBackgroundToTemperature(background(:),ensembleSize,domain(:))
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(background to T) =',ct1-ct0,'sec'
print*,'walltime(background to T) =',wt1-wt0,'sec'


print*,'Converting background to sounding...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call convertBackgroundToSounding(background(:),ensembleSize,domain(:),sounding)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(H of sounding) =',ct1-ct0,'sec'
print*,'walltime(H of sounding) =',wt1-wt0,'sec'


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

!
!  Assimilation on mass grid.
!
print*,'Mapping observations to each mass grid...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call mapObsToEachMassGrid(obsListOfEachMassGrid,sounding,domain_mean)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time (mapObsToEachMassGrid)=',ct1-ct0,'sec'
print*,'walltime (mapObsToEachMassGrid)=',wt1-wt0,'sec'
print*,'There are',count(obsListOfEachMassGrid(:,:,:)%vectorSize>0),'grids may need to be updated.'
print*,'Max obs per grid=',maxval(obsListOfEachMassGrid(:,:,:)%vectorSize)
print*,'Total obs for all grids=',sum(obsListOfEachMassGrid(:,:,:)%vectorSize)


print*,'Starting assimilation...'
allocate( analysis(ensembleSize) )
wt0 = omp_get_wtime()
call cpu_time(ct0)
call assimilate_massGrid(background(:),analysis(:),ensembleSize,domain(:),domain_mean,sounding,obsListOfEachMassGrid)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(assimilation on mass grid) =',ct1-ct0,'sec'
print*,'walltime(assimilation on mass grid) =',wt1-wt0,'sec'

!
!  Assimilation on u grid.
!
print*,'Mapping observations to each u grid...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call mapObsToEachUGrid(obsListOfEachUGrid,sounding,domain_mean)
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
call assimilate_uGrid(background(:),analysis(:),ensembleSize,domain(:),domain_mean,sounding,obsListOfEachUGrid)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(assimilation on u grid) =',ct1-ct0,'sec'
print*,'walltime(assimilation on u grid) =',wt1-wt0,'sec'

!
!  Assimilation on v grid.
!
print*,'Mapping observations to each v grid...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call mapObsToEachVGrid(obsListOfEachVGrid,sounding,domain_mean)
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
call assimilate_vGrid(background(:),analysis(:),ensembleSize,domain(:),domain_mean,sounding,obsListOfEachVGrid)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(assimilation on v grid) =',ct1-ct0,'sec'
print*,'walltime(assimilation on v grid) =',wt1-wt0,'sec'

!
!  Assimilation on w grid.
!
print*,'Mapping observations to each w grid...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call mapObsToEachWGrid(obsListOfEachWGrid,sounding,domain_mean)
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
call assimilate_wGrid(background(:),analysis(:),ensembleSize,domain(:),domain_mean,sounding,obsListOfEachWGrid)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(assimilation on w grid) =',ct1-ct0,'sec'
print*,'walltime(assimilation on w grid) =',wt1-wt0,'sec'

!iz=availableFileID()
!open(iz,file='obsListOfEachUGrid')
!write(iz,*) domain(1)%lon_u(90+10,75+10),domain(1)%lat_u(90+10,75+10)
!do io = 1,obsListOfEachUGrid(90+10,75+10,10)%vectorSize
!    write(iz,*) sounding%obs(obsListOfEachUGrid(90+10,75+10,10)%vector(io))%lon,sounding%obs(obsListOfEachUGrid(90+10,75+10,10)%vector(io))%lat
!enddo
!stop

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
