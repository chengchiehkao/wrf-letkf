
module derivedType


implicit none

type domainInfo
    real(kind=8),pointer,dimension(:,:)   :: lon=>null() , lat=>null()
    real(kind=8),pointer,dimension(:,:)   :: lon_u=>null() , lat_u=>null()
    real(kind=8),pointer,dimension(:,:)   :: lon_v=>null() , lat_v=>null()
    real(kind=8),pointer,dimension(:,:,:) :: pressure=>null() , GPH=>null() , pressure_u=>null() , pressure_v=>null() , pressure_w=>null()
    integer                               :: size_westToEast      , size_southToNorth      , size_bottomToTop
    integer                               :: size_westToEast_stag , size_southToNorth_stag , size_bottomToTop_stag
end type domainInfo

type systemParameter
    integer :: ensembleSize
    logical :: use_sound , use_airep , use_synop , use_amv , use_gpsro , use_airs , use_quikscat , use_ascat , use_iasi , use_oscat
    integer :: varListSize_sound , varListSize_airep , varListSize_synop , varListSize_amv , varListSize_gpsro , varListSize_airs ,varListSize_quikscat , varListSize_ascat , varListSize_iasi , varListSize_oscat
    character(len=10),pointer,dimension(:) :: varList_sound=>null() , varList_airep=>null() , varList_synop=>null() , varList_amv=>null() , varList_gpsro=>null() , varList_airs=>null() ,   &
                                              varList_quikscat=>null() , varList_ascat=>null() , varList_iasi=>null() , varList_oscat=>null()
    logical,pointer,dimension(:) :: use_varList_sound=>null() , use_varList_airep=>null() , use_varList_synop=>null() , use_varList_amv=>null() , use_varList_gpsro=>null() , use_varList_airs=>null() ,  &
                                     use_varList_quikscat=>null() , use_varList_ascat=>null() , use_varList_iasi=>null() , use_varList_oscat=>null()
    real(kind=8) :: rd   , rc    ! Decorrelated & highly correlated distance on horizontal space.
    real(kind=8) :: rd_z , rc_z  ! Decorrelated & highly correlated distance on vertical space.
    real(kind=8) :: inflationFactor
end type systemParameter

type obsChild
    real(kind=8)      :: lon , lat , z , value , error
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
    real(kind=8),pointer,dimension(:,:)   :: mu=>null() , u10=>null() , v10=>null()
    real(kind=8),pointer,dimension(:,:,:) :: u=>null() , v=>null() , w=>null() , t=>null() , qvapor=>null()
    real(kind=8),pointer,dimension(:,:,:) :: ph=>null()  ! PH as perturbation geopotential in WRF.
    real(kind=8),pointer,dimension(:,:,:) :: normalT=>null() , stratifiedMU=>null() , RH=>null()
end type backgroundInfo

end module derivedType
