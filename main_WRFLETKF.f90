
include 'mod_derivedType.f90'
!include 'mod_basicUtility.f90'  ! compiled individually for default real as 8-byte.
include 'mod_IOUtility.f90'
include 'mod_systemUtility.f90'
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

type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid=>null()

integer io,iwe,isn,iz

real         ct0,ct1
real(kind=8) wt0,wt1
!================================================

ensembleSize=36

allocate( domain(ensembleSize) )

print*,'Getting Domain...'
call getDomain(domain(:),ensembleSize)
call getMeanDomain(domain(:),ensembleSize,domain_mean)
print*,'Done.'

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

print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavalible.'
print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) unavalible.'
print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) unavalible.'


print*,repeat('=',20)
print*,'Checking if observations inside vertical domain...'
call check_ifObsInsideVerticalDomain(domain(:),ensembleSize,sounding)
call check_ifObsInsideVerticalDomain(domain(:),ensembleSize,amv)
call check_ifObsInsideVerticalDomain(domain(:),ensembleSize,gpsro)

print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavalible.'
print*,'There are ',count(.not.amv%obs(:)%available),'/',amv%obsNum,'amv(s) unavalible.'
print*,'There are ',count(.not.gpsro%obs(:)%available),'/',gpsro%obsNum,'gpsro(s) unavalible.'


call cpu_time(ct1)
print*,'cpu time (from the beginning except getDomain)=',ct1-ct0,'sec'



print*,'Mapping observations to each grid...'
wt0 = omp_get_wtime()
call cpu_time(ct0)
call mapObsToEachGrid(obsListOfEachGrid,sounding,domain_mean)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time (mapObsToEachGrid)=',ct1-ct0,'sec'
print*,'walltime (mapObsToEachGrid)=',wt1-wt0,'sec'


print*,'There are',count(obsListOfEachGrid(:,:,:)%vectorSize>0),'grids may need to be updated.'
print*,'Max obs per grid=',maxval(obsListOfEachGrid(:,:,:)%vectorSize)
print*,'Total obs for all grids=',sum(obsListOfEachGrid(:,:,:)%vectorSize)


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
print*,'There are ',count(.not.sounding%obs(:)%available),'/',sounding%obsNum,'sounding(s) unavalible.'

print*,'Starting assimilation...'
allocate( analysis(ensembleSize) )
wt0 = omp_get_wtime()
call cpu_time(ct0)
call assimilate(background(:),analysis(:),ensembleSize,domain(:),domain_mean,sounding,obsListOfEachGrid)
call cpu_time(ct1)
wt1 = omp_get_wtime()
print*,'Done.'
print*,'cpu time(assimilation) =',ct1-ct0,'sec'
print*,'walltime(assimilation) =',wt1-wt0,'sec'

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
