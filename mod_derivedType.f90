
module derivedType


implicit none

type domainInfo
    real(kind=8),pointer,dimension(:,:)   :: lon=>null() , lat=>null()
    real(kind=8),pointer,dimension(:,:)   :: lon_u=>null() , lat_u=>null()
    real(kind=8),pointer,dimension(:,:)   :: lon_v=>null() , lat_v=>null()
    real(kind=8),pointer,dimension(:,:,:) :: pressure=>null() , GPH=>null() , pressure_w=>null()
    integer                               :: size_westToEast      , size_southToNorth      , size_bottomToTop
    integer                               :: size_westToEast_stag , size_southToNorth_stag , size_bottomToTop_stag
end type domainInfo

type obsChild
    real(kind=8)      :: lon , lat , z , var , error
    real(kind=8)      :: innov
    real(kind=8)      :: rfict  ! for GPSRO
    character(len=10) :: type , instrument , zName , varName
    real(kind=8),pointer,dimension(:) :: background=>null()
    logical           :: available
end type obsChild

type obsParent
    integer :: obsNum = 0
    type(obsChild),pointer,dimension(:) :: obs=>null()
end type obsParent

type integerVector
    integer,pointer,dimension(:) :: vector=>null()
    integer :: vectorSize = 0
end type integerVector

type backgroundInfo
    real(kind=8),pointer,dimension(:,:)   :: mu=>null()
    real(kind=8),pointer,dimension(:,:,:) :: u=>null() , v=>null() , w=>null() , t=>null() , qvapor=>null()
    real(kind=8),pointer,dimension(:,:,:) :: ph=>null()  ! PH as perturbation geopotential in WRF.
    real(kind=8),pointer,dimension(:,:,:) :: normalT=>null() , stratifiedDryAirMass=>null()
end type backgroundInfo

end module derivedType
