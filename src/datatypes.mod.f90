module datatypes
  use, intrinsic :: iso_fortran_env, dp=>real64, sp=>real32, in=>int32
  use md_interface_lm3ppa, only: myinterface
  use md_params_core

! define data types and constants
 implicit none
! ---- public types -------
 public :: spec_data_type, cohort_type, vegn_tile_type
! ------ public subroutines ---------
public :: initialize_PFT_data, initialize_soilpars
public :: Zero_diagnostics, hourly_diagnostics, daily_diagnostics, &
          annual_diagnostics
public :: qscomp, calc_esat


! ---- public variables ---------
 public :: spdata, soilpars
 ! parameters
 public :: MaxCohortID, &
    K1, K2, K_nitrogen, etaN, MLmixRatio, &
    fsc_fine, fsc_wood,  &
    GR_factor,  l_fract, retransN, f_initialBSW, &
    f_N_add, A_mort, B_mort,DBHtp

!===============constants===============
 logical, public, parameter :: read_from_parameter_file = .TRUE.
 integer, public, parameter :: days_per_year  = 365
 integer, public, parameter :: hours_per_year = 365 * 24  ! 8760
 real,    public, parameter :: seconds_per_year = 365. * 24. * 3600.
 real,    public, parameter :: seconds_per_day = 24. * 3600.

 integer, public, parameter :: max_lev  = 3 ! Soil layers, for soil water dynamics
 integer, public, parameter :: num_l    = 3 ! Soil layers
 integer, public, parameter :: LEAF_ON  = 1
 integer, public, parameter :: LEAF_OFF = 0
 integer, public, parameter :: & ! physiology types
 PT_C3        = 0, &
 PT_C4        = 1

 ! Soil water hydrualics
 real, public, parameter :: rzone = 2.0 !m
 real,public, parameter ::  thksl(max_lev)=(/0.05,0.45,1.5/) ! m, thickness of soil layers
 real, public, parameter :: psi_wilt  = -150.0  ! matric head at wilting
 real, public, parameter :: K_rel_min = 1.e-12
 real, public, parameter :: rate_fc   = 0.1/86400 ! 0.1 mm/d drainage rate at FC
 real, public, parameter :: ws0 = 0.02 ! hygroscopic point
 real, public, parameter :: Edepth = 0.05 !m, the depth of soil for surface evaporation
 integer, public, parameter :: & ! soil types
                    Sand        = 1, &
                    LoamySand   = 2, &
                    SandyLoam   = 3, &
                    SiltLoam    = 4, &
                    FrittedClay = 5, &
                    Loam        = 6, &
                    Clay        = 7

 integer, public, parameter :: & ! phenology type
 PHEN_DECIDIOUS = 0, &
 PHEN_EVERGREEN = 1

 ! Soil SOM reference C/N ratios
 real, public, parameter :: CN0metabolicL  = 15.0 ! 25.0 ! 15.0
 real, public, parameter :: CN0structuralL = 40.0 ! 55.0 ! 35.0

! Physical constants
real, public, parameter :: TFREEZE = 273.16
real, public, parameter :: Rugas = 8.314472 ! universal gas constant, J K-1 mol-1
real, public, parameter :: mol_C = 12.0e-3 ! molar mass of carbon, kg
real, public, parameter :: mol_air = 28.96440e-3 ! molar mass of air, kg
real, public, parameter :: mol_CO2 = 44.00995e-3 ! molar mass of CO2,kg
real, public, parameter :: mol_h2o = 18.0e-3 ! molar mass of water, kg
real, public, parameter :: cpair = 1010.
real, public, parameter :: H2OLv0=2.501e6        !latent heat H2O (J/kg)
real, public, parameter :: p_sea = 101325. !1.e5           ! atmospheric pressure  (Pa)
real, public, parameter :: DENS_H2O = 1000. ! kg m-3

!===============data types ==============================
!-----------PFT data type----------------
type spec_data_type
  integer :: lifeform     ! 0 for grasses, 1 for trees
  integer :: phenotype    ! phenology type: 0 for deciduous, 1 for evergreen
  integer :: pt           ! photosynthetic physiology of species
  ! leaf traits
  real    :: LMA          ! leaf mass per unit area, kg C/m2
  real    :: leafLS       ! leaf life span
  real    :: alpha_L      ! leaf turn over rate
  real    :: LNA          ! leaf Nitrogen per unit area, kg N/m2
  real    :: LNbase       ! basal leaf Nitrogen per unit area, kg N/m2, (Rubisco)
  real    :: CNleafsupport! leaf structural tissues, 175
  real    :: leaf_size    ! characteristic leaf size
  real    :: alpha_phot   ! photosynthesis efficiency
  real    :: m_cond       ! factor of stomatal conductance
  real    :: Vmax         ! max rubisco rate, mol m-2 s-1
  real    :: Vannual      ! annual productivity per unit area at full fun (kgC m-2 yr-1)
  real    :: gamma_L      ! leaf respiration coeficient (per yr)
  real    :: gamma_LN     ! leaf respiration coeficient per unit N
  real    :: wet_leaf_dreg ! wet leaf photosynthesis down-regulation
  ! root traits
  real    :: rho_FR       ! material density of fine roots (kgC m-3)
  real    :: root_r       ! radius of the fine roots, m
  real    :: root_zeta    ! e-folding parameter of root vertical distribution (m)
  real    :: root_frac(max_lev)    ! root fraction
  real    :: SRA          ! speific fine root area, m2/kg C
  real    :: SRL          ! specific root lenght
  real    :: gamma_FR     ! Fine root respiration rate, kgC kgC-1 yr-1
  real    :: alpha_FR     ! Turnover rate of Fine roots, fraction yr-1
  real    :: Kw_root      ! fine root water donductivity mol m m-2 s−1 MPa−1 ! 
  real    :: root_perm
!  real    :: rho_N_up0   ! maximum N uptake rate
!  real    :: N_roots0    ! root biomass at half of max. N-uptake rate
  real    :: NfixRate0    ! Reference N fixation rate (kgN kgC-1 root)
  real    :: NfixCost0    ! Carbon cost of N fixation (kgC kgN-1)
  ! wood traits
  real    :: rho_wood     ! woody density, kg C m-3 wood
  real    :: gamma_SW     ! sapwood respiration rate, kgC m-2 Acambium yr-1
  real    :: taperfactor

  ! Allometry
  real    :: alphaHT, thetaHT ! height = alphaHT * DBH ** thetaHT
  real    :: alphaCA, thetaCA ! crown area = alphaCA * DBH ** thetaCA
  real    :: alphaBM, thetaBM ! biomass = alphaBM * DBH ** thetaBM
  real    :: phiRL            ! ratio of fine root to leaf area
  real    :: phiCSA           ! ratio of sapwood CSA to target leaf area
  real    :: tauNSC           ! residence time of C in NSC (to define storage capacity)
  real    :: fNSNmax          ! multilier for NSNmax

  ! Default C/N ratios
  real    :: CNleaf0
  real    :: CNroot0
  real    :: CNsw0
  real    :: CNwood0
  real    :: CNseed0
  ! phenology
  real    :: tc_crit         ! K, for turning OFF a growth season
  real    :: tc_crit_on      ! K, for turning ON a growth season
  real    :: gdd_crit        ! K, critical value of GDD5 for turning ON growth season
  !  vital rates
  real    :: maturalage       ! the age that can reproduce
  real    :: v_seed           ! fracton of G_SF to G_F
  real    :: seedlingsize     ! size of the seedlings, kgC/indiv
  real    :: prob_g,prob_e    ! germination and establishment probabilities
  real    :: mortrate_d_c     ! yearly mortality rate in canopy
  real    :: mortrate_d_u     ! yearly mortality rate in understory
  ! Population level variables
  real    :: LAImax,underLAImax ! max. LAI
  real    :: LAI_light        ! light controlled maximum LAI
  integer :: n_cc             ! for calculating LAImax via cc%LAImax derived from cc%NSN
  real    :: layerfrac        ! species layer fraction
  real    :: internal_gap_frac ! fraction of internal gaps in the canopy
  ! "internal" gaps are the gaps that are created within the canopy by the
  ! branch fall processes.
end type

!----------cohort type with conserved quantitites-----------------
type :: cohort_conserved_type

end type cohort_conserved_type


!----------cohort type with re-calculated quantitites-----------------
type :: cohort_type

! ----- carbon fluxes
  real :: gpp  = 0.0 ! gross primary productivity kg C/timestep
  real :: npp  = 0.0 ! net primary productivity kg C/timestep
  real :: resp = 0.0 ! plant respiration
  real :: resl = 0.0 ! leaf respiration
  real :: resr = 0.0 ! root respiration
  real :: resg = 0.0 ! growth respiration
  real :: NPPleaf,NPProot,NPPwood ! to record C allocated to leaf, root, and wood
  real :: dailyTrsp
  real :: dailyGPP   ! kgC/tree day-1
  real :: dailyNPP
  real :: dailyResp
  real :: dailyNup
  real :: dailyfixedN
  real :: annualTrsp
  real :: annualGPP ! C flux/tree
  real :: annualNPP
  real :: annualResp

! ---- Nitrogen model related parameters
  real    :: NSNmax = 0.
  real    :: NSN = 0.    ! non-structural N pool
  real    :: leafN = 0.
  real    :: sapwN= 0.
  real    :: woodN = 0. ! N of heart wood
  real    :: rootN = 0. ! N of fine roots
  real    :: seedN = 0. !
  real    :: N_uptake = 0.
  real    :: annualNup  = 0.0
  real    :: fixedN ! fixed N at each stem per tree
  real    :: annualfixedN = 0.0 ! annual N fixation per unit crown area

  ! TODO: see if we can make bl_max, br_max local variables
  real    :: bl_max  = 0.0 ! Max. leaf biomass, kg C/individual
  real    :: br_max  = 0.0 ! Max. fine root biomass, kg C/individual
  real    :: CSAsw   = 0.0
  real    :: topyear = 0.0 ! the years that a plant in top layer
  real    :: DBH_ys        ! DBH at the begining of a year (growing season)

! ---- water uptake-related variables
  real    :: root_length(max_lev) ! m
  real    :: rootarea ! total fine root area per tree
  real    :: rootdepth  ! maximum depth of fine roots
  real    :: rootareaL(max_lev) = 0.0 ! Root length per layer, m of root/m
  real    :: WupL(max_lev) = 0.0 ! normalized vertical distribution of uptake
  real    :: W_supply  ! potential water uptake rate per unit time per tree
  real    :: transp   ! transpiration rate per tree per hour
  real    :: uptake_frac(max_lev) ! for LM3 soil water uptake, Weng, 2017-10-28
  real    :: K_r,r_r
  real    :: root_zeta

! for photosynthesis
  real :: An_op = 0.0 ! mol C/(m2 of leaf per year)
  real :: An_cl = 0.0 ! mol C/(m2 of leaf per year)
  real :: w_scale =-9999
  real :: C_growth = 0.0 ! carbon gain since last growth, kg C/individual
  real :: N_growth = 0.0 ! Nitrogen used for plant tissue growth
  real :: extinct = 0.75     ! light extinction coefficient in the canopy for photosynthesis

end type cohort_type

!-----------conserved variables----------------
! contains quantities that are defined only at tile level (not cohort level), e.g. soil-related variables
type :: vegn_tile_conserved_type
xxx
end type vegn_tile_conserved_type


!-----------re-calculated variables----------------
type :: vegn_tile_type
   integer :: n_cohorts = 0
   integer :: n_years   = 0
   integer :: n_canopycc = 0
   type(cohort_type), pointer :: cohorts(:)=>NULL()
   real :: area  ! m2
   real :: age=0 ! tile age
   real :: nindivs  ! New varaibles at tile level xxx
   real :: DBH ! New varaibles at tile level xxx
   real :: nindivs12
   real :: DBH12
   real :: DBH12pow2
   real :: QMD
   real :: MaxAge
   real :: MaxVolume
   real :: MaxDBH

   ! leaf area index
   real :: LAI  ! leaf area index
   real :: CAI  ! crown area index
   real :: LAIlayer(0:10) = 0.0 ! LAI of each crown layer, max. 9
   ! uptake-related variables
   real :: root_distance(max_lev) ! characteristic half-distance between fine roots, m
   ! averaged quantities for PPA phenology
   real :: tc_daily = 0.0
   real :: gdd      = 0.0 ! growing degree-days
   real :: tc_pheno = 0.0 ! smoothed canopy air temperature for phenology

   ! litter and soil carbon pools
   real :: litter = 0.0 ! litter flux
   real :: MicrobialC  = 0  ! Microbes (kg C/m2)
   real :: metabolicL  = 0  ! fast soil carbon pool, (kg C/m2)
   real :: structuralL = 0  ! slow soil carbon pool, (kg C/m2)

!!  Nitrogen pools, Weng 2014-08-08
   real :: MicrobialN= 0
   real :: metabolicN = 0  ! fast soil nitrogen pool, (kg N/m2)
   real :: structuralN = 0  ! slow soil nitrogen pool, (kg N/m2)
   real :: mineralN = 0.  ! Mineral nitrogen pool, (kg N/m2)
   real :: totN = 0.
   real :: N_input        ! annual N input (kgN m-2 yr-1)
   real :: N_uptake  = 0.  ! kg N m-2 hour-1
   real :: annualN = 0.0  ! annual available N in a year
   real :: Nloss_yr= 0.0  ! annual N loss
   real :: N_P2S_yr= 0.0  ! annual N from plants to soil
   real :: previousN      ! an weighted annual available N
   real :: initialN0

! Soil water
   integer :: soiltype       ! lookup table for soil hydrologic parameters
   real :: FLDCAP,WILTPT  ! soil property: field capacity and wilting point (0.xx)
   real :: evap           ! kg m-2 per unit fast time step (mm/hour)
   real :: transp         ! kg m-2 hour-1
   real :: runoff        ! Water runoff of the veg tile, unit?
   real :: thetaS     ! moisture index (ws - wiltpt)/(fldcap - wiltpt)
   real :: wcl(max_lev)   ! volumetric soil water content for each layer
   real :: soilWater      ! kg m-2 in root zone
! ---- water uptake-related variables
  real    :: RAI ! root area index
  real    :: RAIL(max_lev) = 0.0 ! Root length per layer, m of root/m
  real    :: W_uptake  ! water uptake rate per unit time per m2

!  Carbon fluxes
   real :: gpp =0 ! gross primary production, kgC m-2 yr-1
   real :: npp =0 ! net primary productivity
   real :: resp = 0 ! auto-respiration of plants
   real :: nep =0 ! net ecosystem productivity
   real :: rh  =0 ! soil carbon lost to the atmosphere
   ! daily diagnostics
   real :: dailyGPP
   real :: dailyNPP
   real :: dailyResp
   real :: dailyRh
   real :: dailyNup
   real :: dailyfixedN
   ! for annual diagnostics
   real :: dailyPrcp=0.0, annualPrcp = 0.0 ! mm m-2 yr-1
   real :: dailyTrsp=0.0,dailyEvap=0.0, dailyRoff=0.0 ! mm m-2 yr-1
   real :: annualTrsp=0.0,annualEvap=0.0, annualRoff=0.0 ! mm m-2 yr-1
   real :: annualGPP = 0.0 ! kgC m-2 ground yr-1
   real :: annualNPP = 0.0
   real :: annualResp = 0.0
   real :: annualRh   = 0.0
   real :: annualNup       ! accumulated N uptake kgN m-2 yr-1
   real :: annualfixedN = 0.  ! fixe N in a tile
   ! for annual reporting at tile level
   real :: NSC, SeedC, leafC, rootC, SapwoodC, WoodC
   real :: NSN, SeedN, leafN, rootN, SapwoodN, WoodN
   real :: totSeedC,totSeedN,totNewCC, totNewCN
end type vegn_tile_type

type :: soil_pars_type
  real :: GMD ! geometric mean partice diameter, mm
  real :: GSD ! geometric standard deviation of particle size
  real :: vwc_wilt
  real :: vwc_fc
  real :: vwc_sat
  real :: vlc_min
  real :: k_sat_ref ! hydraulic conductivity of saturated soil, kg/(m2 s)
  real :: psi_sat_ref ! saturation soil water potential, m
  real :: chb         ! Soil texture parameter
  real :: alpha       ! vertical changes of soil property, 1: no change
  real :: heat_capacity_dry
  real::  tfreeze
end type soil_pars_type

type :: soil_prog_type
  real wl
  real ws
  real T
end type soil_prog_type

type :: soil_tile_type
   integer :: tag ! kind of the soil
   type(soil_pars_type) :: pars
   type(soil_prog_type), pointer :: prog(:)
   real,                 pointer :: w_fc(:)
   real,                 pointer :: w_wilt(:)
   real :: Eg_part_ref
   real :: z0_scalar
   ! data that were local to soil.f90
   real, pointer :: heat_capacity_dry(:)
   real, pointer :: e(:),f(:)
   ! added to avoid recalculation of soil hydraulics in case of Darcy uptake
   real          :: uptake_T ! update temperature from previous time step
   real, pointer :: psi(:) ! soil water potential
end type soil_tile_type

! PFT-specific parameters
type(spec_data_type), save :: spdata(0:MSPECIES) ! define PFTs
! Soil
type(soil_pars_type), save :: soilpars(n_dim_soil_types) ! soil parameters

!---------------------------------------
integer :: MaxCohortID = 0

! Constants:
! Soil water properties
! real   :: soiltype = SandyLoam  ! 1 Sand; 2
real   :: FLDCAP = 0.4  ! vol/vol
real   :: WILTPT != 0.05 ! vol/vol
! Carbon pools
real :: K1 != 2 ! Fast soil C decomposition rate (yr-1)
real :: K2 != 0.05 ! slow soil C decomposition rate (yr-1)
real :: K_nitrogen != 8.0     ! mineral Nitrogen turnover rate
real :: etaN       != 0.025   ! N loss through runoff (organic and mineral)
real :: MLmixRatio != 0.8     ! the ratio of C and N returned to litters from microbes
real :: l_fract    != 0.0     ! 0.25  ! 0.5 ! fraction of the carbon retained after leaf drop
real :: retransN   != 0.0     ! retranslocation coefficient of Nitrogen
real :: f_N_add != 0.02       ! re-fill of N for sapwood
real :: f_initialBSW != 0.2   !0.01
real :: LMAmin     = 0.02    ! minimum LMA, boundary condition
real :: fsc_fine   = 1.0     ! fraction of fast turnover carbon in fine biomass
real :: fsc_wood   = 0.0     ! fraction of fast turnover carbon in wood biomass
real :: GR_factor  = 0.33    ! growth respiration factor

! Ensheng's growth parameters:
real :: f_LFR_max =0.85 ! max allocation to leaves and fine roots ! wood_fract_min = 0.15 ! for understory mortality rate is calculated as:
! deathrate = mortrate_d_u * (1+A*exp(B*DBH))/(1+exp(B*DBH))
real :: A_mort     = 9.0    ! A coefficient in understory mortality rate correction, 1/year
real :: B_mort     = -60.0  ! B coefficient in understory mortality rate correction, 1/m
real :: DBHtp      = 2.0    !  m, for canopy tree's mortality rate
! for leaf life span and LMA (leafLS = c_LLS * LMA
real :: c_LLS  = 28.57143   ! yr/ (kg C m-2), 1/LMAs, where LMAs = 0.035

! reduction of bl_max and br_max for the understory vegetation, unitless
real :: understory_lai_factor = 0.25
!real :: rdepth(0: max_lev) = 0.0

! -------- PFT-specific parameters ----------
! c4grass  c3grass  temp-decid  tropical  evergreen  BE  BD  BN  NE  ND  G  D  T  A
! integer :: pt(0:MSPECIES) = 0
!(/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0/) ! 0 for C3, 1 for C4
! integer :: phenotype(0:MSPECIES)= 0
! (/0,  0,  0,  0,  1,  1,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0 /) ! 0 for Deciduous, 1 for evergreen
! integer :: lifeform(0:MSPECIES) = 1 ! life form of PFTs: 0 for grasses, 1 for trees

! root parameters
real :: alpha_FR(0:MSPECIES) = 1.2 ! Fine root turnover rate yr-1
!(/0.8, 0.8,0.8, 0.8, 0.8,0.8,0.8,0.8,1.0,1.0,0.6, 1.0, 0.55, 0.9, 0.55, 0.55/)
real :: rho_FR(0:MSPECIES) = 200 ! woody density, kgC m-3
real :: root_r(0:MSPECIES) = 2.9E-4
!(/1.1e-4, 1.1e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 2.9e-4, 1.1e-4, 1.1e-4, 2.2e-4, 2.2e-4/)
real :: root_zeta(0:MSPECIES) = 0.29 !
real :: Kw_root(0:MSPECIES)= 6.3E-8 * (1000000.0/18.0)*1.e-6 ! mol /(s m2 Mpa) ! 6.3±3.1×10−8 m s−1 MPa−1
!(/1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5/)
   ! fine root membrane permeability per unit membrane area, kg/(m3 s).
   ! Root membrane permeability is "high" for the value from Siqueira et al., 2008,
! Water Resource Research Vol. 44, W01432, converted to mass units
!real :: rho_N_up0(0:MSPECIES) = 0.5 ! fraction of mineral N per hour
!real :: N_roots0(0:MSPECIES) = 0.3 ! kgC m-2

real :: leaf_size(0:MSPECIES)= 0.04 !
! photosynthesis parameters
real :: Vmax(0:MSPECIES)= 35.0E-6 ! mol m-2 s-1
real :: Vannual(0:MSPECIES) = 1.2 ! kgC m-2 yr-1
real :: wet_leaf_dreg(0:MSPECIES) = 0.3 ! wet leaf photosynthesis down-regulation: 0.3 means
        ! photosynthesis of completely wet leaf will be 30% less than that of dry one,
        ! provided everything else is the same
!(/1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2/)
real :: m_cond(0:MSPECIES)= 7.0 !
real :: alpha_phot(0:MSPECIES)=  0.06 !
real :: gamma_L(0:MSPECIES)= 0.02 !
real :: gamma_LN(0:MSPECIES)= 70.5 ! 25.0  ! kgC kgN-1 yr-1
real :: gamma_SW(0:MSPECIES)= 0.08 ! 5.0e-4 ! kgC m-2 Acambium yr-1
real :: gamma_FR(0:MSPECIES)= 12.0 ! 15 !kgC kgN-1 yr-1 ! 0.6: kgC kgN-1 yr-1
real :: tc_crit(0:MSPECIES)= 283.16 ! OFF
real :: tc_crit_on(0:MSPECIES)= 280.16 ! ON
real :: gdd_crit(0:MSPECIES)= 280.0 ! Simulations 280, 240, 200

! Allometry parameters
real :: alphaHT(0:MSPECIES)      = 36.0
real :: thetaHT(0:MSPECIES)      = 0.5 !
real :: alphaCA(0:MSPECIES)      = 150.0
real :: thetaCA(0:MSPECIES)      = 1.5
real :: alphaBM(0:MSPECIES)      = 5200.0
real :: thetaBM(0:MSPECIES)      = 2.5

! Reproduction parameters
real :: maturalage(0:MSPECIES)   = 5.0  ! year
real :: v_seed(0:MSPECIES)       = 0.1  ! fraction of allocation to wood+seeds
real :: seedlingsize(0:MSPECIES) = 0.05 ! kgC
real :: prob_g(0:MSPECIES)       = 1.0
real :: prob_e(0:MSPECIES)       = 1.0

! Mortality
real :: mortrate_d_c(0:MSPECIES) = 0.01 ! yearly
real :: mortrate_d_u(0:MSPECIES) = 0.075

! Leaf parameters
! real :: LMA(0:MSPECIES)         = 0.035  ! (Simulations: 0.035, 0.085, 0.135) leaf mass per unit area, kg C/m2
!(/0.04,    0.04,    0.035,   0.035,   0.140,  0.032, 0.032,  0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036,   0.036  /)
real :: leafLS(0:MSPECIES) = 1.0
real :: LNbase(0:MSPECIES)        = 0.8E-3 !functional nitrogen per unit leaf area, kg N/m2
real :: CNleafsupport(0:MSPECIES) = 80.0 ! CN ratio of leaf supporting tissues
real :: rho_wood(0:MSPECIES)      = 300.0 ! kgC m-3 (Simulations: 300, 600, 800)
real :: taperfactor(0:MSPECIES)   = 0.75 ! taper factor, from a cylinder to a tree
real :: LAImax(0:MSPECIES)        = 3.5 ! maximum LAI for a tree
real :: LAI_light(0:MSPECIES)     = 4.0 ! maximum LAI limited by light
real :: tauNSC(0:MSPECIES)        = 3 ! 3 ! NSC residence time,years
real :: fNSNmax(0:MSPECIES)       = 5 ! 5 ! multilier for NSNmax as sum of potential bl and br
! real :: phiRL(0:MSPECIES)       = 3.5 ! ratio of fine root area to leaf area (Root:Shoot ratio simulations: 3.5, 5, 7)
real :: phiCSA(0:MSPECIES)        = 0.25E-4 ! ratio of sapwood area to leaf area
! C/N ratios for plant pools
real :: CNleaf0(0:MSPECIES)   = 25. ! C/N ratios for leaves
real :: CNsw0(0:MSPECIES)     = 350.0 ! C/N ratios for woody biomass
real :: CNwood0(0:MSPECIES)   = 350.0 ! C/N ratios for woody biomass
real :: CNroot0(0:MSPECIES)   = 40.0 ! C/N ratios for leaves ! Gordon & Jackson 2000
real :: CNseed0(0:MSPECIES)   = 20.0 ! C/N ratios for seeds
real :: NfixRate0(0:MSPECIES) = 0.0 !Reference N fixation rate (0.03 kgN kgC-1 root yr-1)
real :: NfixCost0(0:MSPECIES) = 12.0 ! FUN model, Fisher et al. 2010, GBC
real :: internal_gap_frac(0:MSPECIES)= 0.1 ! The gaps between trees

! soil parameters (passed through R)
! Coarse  Medium   Fine    CM     CF     MF    CMF    Peat    MCM
  ! real :: GMD(n_dim_soil_types) = & ! geometric mean partice diameter, mm
  ! (/ 0.7, 0.4, 0.3, 0.1, 0.1, 0.07, 0.007, 0.3, 0.3 /)
  ! real :: GSD(n_dim_soil_types) = & ! geometric standard deviation of particle size
  ! (/5.0, 5.3, 7.4, 6.1, 6.1, 14.0, 15.0, 7.4, 7.4 /)
  ! real :: vwc_sat(n_dim_soil_types)= &
  !  (/ 0.380, 0.445, 0.448, 0.412, 0.414, 0.446, 0.424, 0.445, 0.445   /)
  ! !real :: vlc_min(n_dim_soil_types)
  ! real :: k_sat_ref(n_dim_soil_types)= & ! mol/(s MPa m) , hydraulic conductivity of saturated soil,
  ! (/ 130.8, 75.1, 53.2, 12.1, 11.1, 12.7, 1.69, 53.2, 53.2 /)
  ! real :: psi_sat_ref(n_dim_soil_types) = & ! Pa
  ! (/ -600., -790., -910., -1580., -1680., -1880., -5980., -790., -790./)
  ! real :: chb(n_dim_soil_types) = &         ! Soil texture parameter
  ! (/   3.5,   6.4,  11.0,   4.8,   6.3,   8.4,   6.3,   6.4,   6.4   /)
  ! real :: alphaSoil(n_dim_soil_types) = 1.0       ! *** REPLACE LATER BY alpha(layer)
  ! real :: heat_capacity_dry(n_dim_soil_types) = &
  ! (/ 1.2e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.1e6, 1.4e6,   1.0   /)

!----- Initial conditions -------------
integer :: init_n_cohorts                        = MAX_INIT_COHORTS
! integer :: init_cohort_species(MAX_INIT_COHORTS) = 2
! real    :: init_cohort_nindivs(MAX_INIT_COHORTS) = 1.0  ! initial individual density, individual/m2
real    :: init_cohort_bl(MAX_INIT_COHORTS)      = 0.0  ! initial biomass of leaves, kg C/individual
real    :: init_cohort_br(MAX_INIT_COHORTS)      = 0.0  ! initial biomass of fine roots, kg C/individual
! real    :: init_cohort_bsw(MAX_INIT_COHORTS)     = 0.05 ! initial biomass of sapwood, kg C/individual
! real    :: init_cohort_bHW(MAX_INIT_COHORTS)     = 0.0  ! initial biomass of heartwood, kg C/tree
real    :: init_cohort_seedC(MAX_INIT_COHORTS)   = 0.0  ! initial biomass of seeds, kg C/individual
! real    :: init_cohort_nsc(MAX_INIT_COHORTS)     = 0.05 ! initial non-structural biomass, kg C/
!  initial soil Carbon and Nitrogen for a vegn tile, Weng 2012-10-24 (pass them through R)
! (passed through R)
! real   :: init_fast_soil_C  = 0.0  ! initial fast soil C, kg C/m2
! real   :: init_slow_soil_C  = 0.0  ! initial slow soil C, kg C/m2
! real   :: init_Nmineral = 0.015  ! Mineral nitrogen pool, (kg N/m2)
real   :: N_input   ! = 0.0008 ! annual N input to soil N pool, kgN m-2 yr-1

! (passed through R)
! logical   :: outputhourly = .False.
! logical   :: outputdaily  = .True.
! logical   :: do_U_shaped_mortality = .False.
logical   :: update_annualLAImax != .False.
! logical   :: do_closedN_run = .True. !.False.

 contains
!=============== subroutines =================================

! ================Parameter initialization ===================
! =========================================================================
subroutine initialize_soilpars(namelistfile)
  use md_interface_lm3ppa, only: myinterface
   character(len=50),intent(in) :: namelistfile
  
  ! initialize soil parameters
    soilpars%GMD               = myinterface%params_soil%GMD(:) ! geometric mean partice diameter, mm
    soilpars%GSD               = myinterface%params_soil%GSD(:) ! geometric standard deviation of particle size
    soilpars%vwc_sat           = myinterface%params_soil%vwc_sat(:)
    soilpars%k_sat_ref         = myinterface%params_soil%k_sat_ref(:) ! hydraulic conductivity of saturated soil, kg/(m2 s)
    soilpars%psi_sat_ref       = myinterface%params_soil%psi_sat_ref(:) ! saturation soil water potential, m
    soilpars%chb               = myinterface%params_soil%chb(:)       ! Soil texture parameter
    soilpars%alpha             = myinterface%params_soil%alphaSoil(:)       ! *** REPLACE LATER BY alpha(layer)
    soilpars%heat_capacity_dry = myinterface%params_soil%heat_capacity_dry(:)

  ! ---- derived constant soil parameters
  ! w_fc (field capacity) set to w at which hydraulic conductivity equals
  ! a nominal drainage rate "rate_fc"
  ! w_wilt set to w at which psi is psi_wilt
  soilpars%vwc_wilt = soilpars%vwc_sat &
          *(soilpars%psi_sat_ref/(psi_wilt*soilpars%alpha))**(1/soilpars%chb)
  soilpars%vwc_fc = soilpars%vwc_sat &
              *(rate_fc/(soilpars%k_sat_ref*soilpars%alpha**2))**(1/(3+2*soilpars%chb))
  soilpars%vlc_min = soilpars%vwc_sat*K_rel_min**(1/(3+2*soilpars%chb))

end subroutine initialize_soilpars

! ================================================
subroutine initialize_PFT_data() !namelistfile

  use md_interface_lm3ppa, only: myinterface

! Initialize PFT parameters
   ! character(len=50),intent(in) :: namelistfile

  ! ---- local vars
  integer :: i

! !  Read parameters from the parameter file (namelist)
!   if(read_from_parameter_file)then
!       nml_unit = 999
!       open(nml_unit, file=namelistfile, form='formatted', action='read', status='old')
!       read (nml_unit, nml=vegn_parameters_nml, iostat=io, end=10)
! 10    close (nml_unit)
!    endif
!       write(*,nml=vegn_parameters_nml)

  ! initialize vegetation data structure
  spdata%pt            = myinterface%params_species(:)%pt
  spdata%phenotype     = myinterface%params_species(:)%phenotype !phenotype  
  spdata%Vmax          = Vmax
  spdata%Vannual       = Vannual
  spdata%m_cond        = m_cond
  spdata%alpha_phot    = alpha_phot
  spdata%wet_leaf_dreg = wet_leaf_dreg
  spdata%gamma_L       = gamma_L
  spdata%gamma_LN      = gamma_LN
  spdata%gamma_SW      = gamma_SW
  spdata%gamma_FR      = gamma_FR

  spdata%rho_FR    = rho_FR
  spdata%root_r    = root_r
  spdata%root_zeta = root_zeta
  spdata%Kw_root   = Kw_root
!  spdata%rho_N_up0 = rho_N_up0
!  spdata%N_roots0  = N_roots0

  spdata%leaf_size = leaf_size
  spdata%tc_crit   = tc_crit
  spdata%gdd_crit  = gdd_crit

! Plant traits
  spdata%LMA           = myinterface%params_species(:)%LMA ! leaf mass per unit area, kg C/m2
  spdata%LNbase        = myinterface%params_species(:)%LNbase   ! Basal leaf nitrogen per unit area, kg N/m2
  spdata%CNleafsupport = CNleafsupport
  spdata%lifeform      = myinterface%params_species(:)%lifeform
  spdata%alphaHT       = alphaHT
  spdata%thetaHT       = thetaHT
  spdata%alphaCA       = alphaCA
  spdata%thetaCA       = thetaCA
  spdata%alphaBM       = alphaBM
  spdata%thetaBM       = thetaBM

  spdata%maturalage   = myinterface%params_species(:)%maturalage
  spdata%v_seed       = v_seed
  spdata%seedlingsize = myinterface%params_species(:)%seedlingsize
  spdata%prob_g       = prob_g
  spdata%prob_e       = prob_e
  spdata%mortrate_d_c = myinterface%params_species(:)%mortrate_d_c
  spdata%mortrate_d_u = myinterface%params_species(:)%mortrate_d_u
  spdata%rho_wood     = rho_wood
  spdata%taperfactor  = taperfactor
  spdata%laimax       = myinterface%params_species(:)%laimax
  spdata%underLAImax  = laimax
  spdata%LAI_light    = myinterface%params_species(:)%LAI_light
  spdata%tauNSC       = tauNSC
  spdata%fNSNmax      = fNSNmax
  spdata%phiRL        = myinterface%params_species(:)%phiRL
  spdata%phiCSA       = myinterface%params_species(:)%phiCSA
  ! root urnover rate
  spdata%alpha_FR     = alpha_FR


!! Nitrogen Weng 2012-10-24
! spdata%CNleaf0 = CNleaf0
  spdata%CNsw0     = CNsw0
  spdata%CNwood0   = CNwood0
  spdata%CNroot0   = CNroot0
  spdata%CNseed0   = CNseed0
  spdata%Nfixrate0 = myinterface%params_species(:)%Nfixrate0
  spdata%NfixCost0 = myinterface%params_species(:)%NfixCost0

  spdata%internal_gap_frac = internal_gap_frac
  do i = 0, MSPECIES
     call init_derived_species_data(spdata(i))
  enddo
  end subroutine initialize_pft_data

!------------------------------------------
 subroutine init_derived_species_data(sp)
   type(spec_data_type), intent(inout) :: sp
   ! ---- local vars ------
   integer :: j
   real :: rdepth(0:max_lev)
   real :: residual

   ! specific root area
   sp%SRA = 2.0/(sp%root_r*sp%rho_FR) ! m2/kgC
   ! root vertical profile
   rdepth=0.0
   do j=1,max_lev
     rdepth(j) = rdepth(j-1)+thksl(j)
     sp%root_frac(j) = exp(-rdepth(j-1)/sp%root_zeta)- &
                       exp(-rdepth(j)  /sp%root_zeta)
   enddo
   residual = exp(-rdepth(max_lev)/sp%root_zeta)
   do j=1,max_lev
      sp%root_frac(j) = sp%root_frac(j) + residual*thksl(j)/rdepth(max_lev)
   enddo

   ! calculate alphaBM parameter of allometry. note that rho_wood was re-introduced for this calculation
   sp%alphaBM    = sp%rho_wood * sp%taperfactor * PI/4. * sp%alphaHT ! 5200
   print*, sp%rho_wood

!  Vmax as a function of LNbase
   sp%Vmax = 0.02 * sp%LNbase ! 0.03125 * sp%LNbase ! Vmax/LNbase= 25E-6/0.8E-3 = 0.03125 !
!  CN0 of leaves
   sp%LNA     = sp%LNbase +  sp%LMA/sp%CNleafsupport
   sp%CNleaf0 = sp%LMA/sp%LNA
!  Leaf life span as a function of LMA
   sp%leafLS = c_LLS * sp%LMA
   if(sp%leafLS>1.0)then
      sp%phenotype = 1
   else
      sp%phenotype = 0
   endif
!  Leaf turnover rate, (leaf longevity as a function of LMA)
   sp%alpha_L = 1.0/sp%leafLS * sp%phenotype

 end subroutine init_derived_species_data

! ============================================================
subroutine qscomp(T, p, qsat)
  real, intent(in) :: T    ! temperature, degK
  real, intent(in) :: p    ! pressure, Pa
  real, intent(out):: qsat ! saturated specific humidity, kg/kg
  !--------local var
  real :: myesat ! sat. water vapor pressure
  real :: Temp ! degC

  ! calculate saturated specific humidity
  Temp = T - 273.16 ! degC
  myesat=MIN(610.78*exp(17.27*Temp/(Temp+237.3)), p) ! Pa

  qsat = 0.622*myesat /(p - 0.378*myesat )

end subroutine qscomp

FUNCTION calc_esat(T) result( out_esat ) ! pressure, Pa
   IMPLICIT NONE
   REAL :: out_esat
   REAL, INTENT(IN) :: T ! degC
   out_esat=610.78*exp(17.27*T/(T+237.3))
END FUNCTION calc_esat
! ==================================

!==============for diagnostics============================================
! Weng, 2016-11-28
subroutine Zero_diagnostics(vegn)
! for annual update
  type(vegn_tile_type), intent(inout) :: vegn
  !-------local var
  type(cohort_type),pointer :: cc
  integer :: i
  !daily
  vegn%dailyfixedN = 0.
  vegn%dailyPrcp   = 0.0
  vegn%dailyTrsp   = 0.0
  vegn%dailyEvap   = 0.0
  vegn%dailyRoff   = 0.0
  vegn%dailyNup    = 0.0
  vegn%dailyGPP    = 0.0
  vegn%dailyNPP    = 0.0
  vegn%dailyResp   = 0.0
  vegn%dailyRh     = 0.0

  !annual
  vegn%annualfixedN = 0.
  vegn%annualPrcp   = 0.0
  vegn%annualTrsp   = 0.0
  vegn%annualEvap   = 0.0
  vegn%annualRoff   = 0.0
  vegn%annualGPP    = 0.0
  vegn%annualNPP    = 0.0
  vegn%annualResp   = 0.0
  vegn%annualRh     = 0.0
  vegn%N_P2S_yr     = 0.
  vegn%annualN      = 0.
  vegn%Nloss_yr     = 0.
  vegn%annualNup    = 0.0

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     cc%C_growth     = 0.0
     cc%N_growth     = 0.0
     cc%gpp          = 0.0
     cc%npp          = 0.0
     cc%resp         = 0.0
     cc%resl         = 0.0
     cc%resr         = 0.0
     cc%resg         = 0.0
     cc%transp       = 0.0
     !daily
     cc%dailyTrsp    = 0.0
     cc%dailyGPP     = 0.0
     cc%dailyNPP     = 0.0
     cc%dailyResp    = 0.0
     cc%dailyNup     = 0.0
     cc%dailyfixedN  = 0.0
     ! annual
     cc%annualTrsp   = 0.0
     cc%annualGPP    = 0.0
     cc%annualNPP    = 0.0
     cc%annualResp   = 0.0
     cc%annualNup    = 0.0
     cc%annualfixedN = 0.0
     cc%NPPleaf      = 0.0
     cc%NPProot      = 0.0
     cc%NPPwood      = 0.0
     cc%DBH_ys       = cc%DBH
  enddo
  
end subroutine Zero_diagnostics

! ========================
subroutine summarize_tile(vegn)
! for annual update
  type(vegn_tile_type), intent(inout) :: vegn
  !-------local var
  type(cohort_type),pointer :: cc
  integer :: i

  ! State variables
  vegn%NSC     = 0.0
  vegn%SeedC   = 0.0
  vegn%leafC   = 0.0
  vegn%rootC   = 0.0
  vegn%SapwoodC= 0.0
  vegn%WoodC   = 0.0

  vegn%NSN     = 0.0
  vegn%SeedN   = 0.0
  vegn%leafN   = 0.0
  vegn%rootN   = 0.0
  vegn%SapwoodN= 0.0
  vegn%WoodN   = 0.0

  vegn%LAI     = 0.0
  vegn%CAI     = 0.0

  ! New tile outputs xxx
  vegn%nindivs    = 0.0
  vegn%DBH        = 0.0
  vegn%nindivs12  = 0.0
  vegn%DBH12      = 0.0
  vegn%DBH12pow2  = 0.0
  vegn%MaxAge     = 0.0
  vegn%MaxVolume  = 0.0
  vegn%MaxDBH     = 0.0

  do i = 1, vegn%n_cohorts
        cc => vegn%cohorts(i)

        ! Vegn C pools:
        vegn%NSC     = vegn%NSC      + cc%NSC       * cc%nindivs
        vegn%SeedC   = vegn%SeedC    + cc%seedC     * cc%nindivs
        vegn%leafC   = vegn%leafC    + cc%bl        * cc%nindivs
        vegn%rootC   = vegn%rootC    + cc%br        * cc%nindivs
        vegn%SapwoodC= vegn%SapwoodC + cc%bsw       * cc%nindivs
        vegn%woodC   = vegn%woodC    + cc%bHW       * cc%nindivs
        vegn%CAI     = vegn%CAI      + cc%crownarea * cc%nindivs
        vegn%LAI     = vegn%LAI      + cc%leafarea  * cc%nindivs
        ! Vegn N pools
        vegn%NSN     = vegn%NSN      + cc%NSN       * cc%nindivs
        vegn%SeedN   = vegn%SeedN    + cc%seedN     * cc%nindivs
        vegn%leafN   = vegn%leafN    + cc%leafN     * cc%nindivs
        vegn%rootN   = vegn%rootN    + cc%rootN     * cc%nindivs
        vegn%SapwoodN= vegn%SapwoodN + cc%sapwN     * cc%nindivs
        vegn%woodN   = vegn%woodN    + cc%woodN     * cc%nindivs
        ! New tile outputs xxx
        vegn%DBH     = vegn%DBH      + cc%dbh       * cc%nindivs
        vegn%nindivs = vegn%nindivs  + cc%nindivs

        if (cc%age    > vegn%MaxAge)       vegn%MaxAge    = cc%age
        if (cc%Volume > vegn%MaxVolume)    vegn%MaxVolume = cc%Volume
        if (cc%DBH    > vegn%MaxDBH)       vegn%MaxDBH    = cc%DBH

        if (cc%dbh > 0.12) then
        vegn%DBH12      = vegn%DBH12     + cc%dbh      * cc%nindivs 
        vegn%nindivs12  = vegn%nindivs12 + cc%nindivs
        vegn%DBH12pow2  = vegn%DBH12pow2 + cc%dbh      * cc%dbh     * cc%nindivs
        
        endif

  enddo

    if (vegn%nindivs>0.0) vegn%DBH     = vegn%DBH / vegn%nindivs  
    if (vegn%nindivs12>0.0) vegn%DBH12 = vegn%DBH12 / vegn%nindivs12  ! vegn%nindivs12 could be zero if all dbh<0.12
    if (vegn%nindivs12>0.0) vegn%QMD   = sqrt(vegn%DBH12pow2 / vegn%nindivs12)  

end subroutine summarize_tile

!=========================================================================
! Hourly fluxes sum to daily
  subroutine hourly_diagnostics(vegn, forcing, iyears, idoy, ihour, out_hourly_tile)

  use md_forcing_lm3ppa, only: climate_type
  use md_interface_lm3ppa, only: outtype_hourly_tile, myinterface

  type(vegn_tile_type), intent(inout) :: vegn
  type(climate_type),intent(in):: forcing
  integer, intent(in) :: iyears, idoy, ihour
  type(outtype_hourly_tile),intent(out) :: out_hourly_tile !!

  !-------local var ------
  type(cohort_type), pointer :: cc    ! current cohort
  integer :: i
  integer :: ntstepsyear !differ

  vegn%age = vegn%age + myinterface%dt_fast_yr
  ! Tile summary
  vegn%GPP    = 0.
  vegn%NPP    = 0.; vegn%Resp   = 0.

  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! cohort daily
     cc%dailyTrsp = cc%dailyTrsp + cc%transp ! kg day-1
     cc%dailyGPP  = cc%dailygpp  + cc%gpp ! kg day-1
     cc%dailyNPP  = cc%dailyNpp  + cc%Npp ! kg day-1
     cc%dailyResp = cc%dailyResp + cc%Resp ! kg day-1

     ! Tile hourly
     vegn%GPP    = vegn%GPP    + cc%gpp    * cc%nindivs
     vegn%NPP    = vegn%NPP    + cc%Npp    * cc%nindivs
     vegn%Resp   = vegn%Resp   + cc%Resp   * cc%nindivs
  enddo
  ! NEP is equal to NNP minus soil respiration
  vegn%nep = vegn%npp - vegn%rh ! kgC m-2 hour-1; time step is hourly
  !! Output horly diagnostics

  ! xxx test: removing condition
    if (.not. myinterface%steering%spinup) then

        out_hourly_tile%year      =  iyears    
        out_hourly_tile%doy       =  idoy   
        out_hourly_tile%hour      =  ihour    
        out_hourly_tile%rad       =  forcing%radiation    !forcingData 
        out_hourly_tile%Tair      =  forcing%Tair         !forcingData  
        out_hourly_tile%Prcp      =  forcing%rain         !forcingData 
        out_hourly_tile%GPP       =  vegn%GPP  
        out_hourly_tile%Resp      =  vegn%resp   
        out_hourly_tile%Transp    =  vegn%transp
        out_hourly_tile%Evap      =  vegn%evap   
        out_hourly_tile%Runoff    =  vegn%runoff   
        out_hourly_tile%Soilwater =  vegn%soilwater
        out_hourly_tile%wcl       =  vegn%wcl(1)    
        out_hourly_tile%FLDCAP    =  vegn%FLDCAP
        out_hourly_tile%WILTPT    =  vegn%WILTPT

    end if

  ! Daily summary:
  vegn%dailyNup  = vegn%dailyNup  + vegn%N_uptake
  vegn%dailyGPP  = vegn%dailyGPP  + vegn%gpp
  vegn%dailyNPP  = vegn%dailyNPP  + vegn%npp
  vegn%dailyResp = vegn%dailyResp + vegn%resp
  vegn%dailyRh   = vegn%dailyRh   + vegn%rh
  vegn%dailyTrsp = vegn%dailyTrsp + vegn%transp
  vegn%dailyEvap = vegn%dailyEvap + vegn%evap
  vegn%dailyRoff = vegn%dailyRoff + vegn%runoff
  vegn%dailyPrcp = vegn%dailyPrcp + forcing%rain * myinterface%step_seconds

  ! print*,'hourly_diagnostics() : vegn%evap, vegn%dailyEvap ', vegn%evap, vegn%dailyEvap

end subroutine hourly_diagnostics

!============================================
  subroutine daily_diagnostics(vegn, iyears, idoy, out_daily_cohorts, out_daily_tile)

    use md_forcing_lm3ppa, only: climate_type
    use md_interface_lm3ppa, only: outtype_daily_cohorts, outtype_daily_tile

    type(vegn_tile_type), intent(inout) :: vegn
    integer, intent(in) :: iyears, idoy
    type(outtype_daily_cohorts), dimension(out_max_cohorts), intent(out) :: out_daily_cohorts
    type(outtype_daily_tile), intent(out) :: out_daily_tile

    !-------local var ------
    type(cohort_type), pointer :: cc    ! current cohort
    integer :: i

    ! re-initialise to avoid elements not updated when number 
    ! of cohorts declines from one year to the next

    if (.not. myinterface%steering%spinup) then 
      out_daily_cohorts(:)%year    = dummy
      out_daily_cohorts(:)%doy     = dummy
      out_daily_cohorts(:)%hour    = dummy
      out_daily_cohorts(:)%cID     = dummy
      out_daily_cohorts(:)%PFT     = dummy
      out_daily_cohorts(:)%layer   = dummy
      out_daily_cohorts(:)%density = dummy
      out_daily_cohorts(:)%f_layer = dummy
      out_daily_cohorts(:)%LAI     = dummy
      out_daily_cohorts(:)%gpp     = dummy
      out_daily_cohorts(:)%resp    = dummy
      out_daily_cohorts(:)%transp  = dummy
      out_daily_cohorts(:)%NPPleaf = dummy
      out_daily_cohorts(:)%NPProot = dummy
      out_daily_cohorts(:)%NPPwood = dummy
      out_daily_cohorts(:)%NSC     = dummy
      out_daily_cohorts(:)%seedC   = dummy
      out_daily_cohorts(:)%leafC   = dummy
      out_daily_cohorts(:)%rootC   = dummy
      out_daily_cohorts(:)%SW_C    = dummy
      out_daily_cohorts(:)%HW_C    = dummy
      out_daily_cohorts(:)%NSN     = dummy
      out_daily_cohorts(:)%seedN   = dummy
      out_daily_cohorts(:)%leafN   = dummy
      out_daily_cohorts(:)%rootN   = dummy
      out_daily_cohorts(:)%SW_N    = dummy
      out_daily_cohorts(:)%HW_N    = dummy
    endif


    ! cohorts output
    do i = 1, vegn%n_cohorts

      cc => vegn%cohorts(i)
      
      if (.not. myinterface%steering%spinup) then 
        out_daily_cohorts(i)%year    = iyears
        out_daily_cohorts(i)%doy     = idoy
        out_daily_cohorts(i)%hour    = i !1.0
        out_daily_cohorts(i)%cID     = cc%ccID
        out_daily_cohorts(i)%PFT     = cc%species
        out_daily_cohorts(i)%layer   = cc%layer
        out_daily_cohorts(i)%density = cc%nindivs*10000
        out_daily_cohorts(i)%f_layer = cc%layerfrac
        out_daily_cohorts(i)%LAI     = cc%LAI
        out_daily_cohorts(i)%gpp     = cc%dailygpp
        out_daily_cohorts(i)%resp    = cc%dailyresp
        out_daily_cohorts(i)%transp  = cc%dailytrsp
        out_daily_cohorts(i)%NPPleaf = cc%NPPleaf
        out_daily_cohorts(i)%NPProot = cc%NPProot
        out_daily_cohorts(i)%NPPwood = cc%NPPwood
        out_daily_cohorts(i)%NSC     = cc%NSC
        out_daily_cohorts(i)%seedC   = cc%seedC
        out_daily_cohorts(i)%leafC   = cc%bl
        out_daily_cohorts(i)%rootC   = cc%br
        out_daily_cohorts(i)%SW_C    = cc%bsw
        out_daily_cohorts(i)%HW_C    = cc%bHW
        out_daily_cohorts(i)%NSN     = cc%NSN*1000
        out_daily_cohorts(i)%seedN   = cc%seedN*1000
        out_daily_cohorts(i)%leafN   = cc%leafN*1000
        out_daily_cohorts(i)%rootN   = cc%rootN*1000
        out_daily_cohorts(i)%SW_N    = cc%sapwN*1000
        out_daily_cohorts(i)%HW_N    = cc%woodN*1000
      endif

      ! annual sum
      cc%annualGPP = cc%annualGPP + cc%dailyGPP
      cc%annualNPP = cc%annualNPP + cc%dailyNPP
      cc%annualResp = cc%annualResp + cc%dailyResp
      cc%annualTrsp = cc%annualTrsp + cc%dailyTrsp
      ! Zero Daily variables
      cc%dailyTrsp = 0.0
      cc%dailyGPP = 0.0
      cc%dailyNPP = 0.0
      cc%dailyResp = 0.0

    enddo

    ! Tile level, daily
    call summarize_tile(vegn)

    if (.not. myinterface%steering%spinup) then 

      out_daily_tile%year      = iyears
      out_daily_tile%doy       = idoy
      out_daily_tile%Tc        = vegn%tc_daily
      out_daily_tile%Prcp      = vegn%dailyPrcp
      out_daily_tile%totWs     = vegn%soilwater
      out_daily_tile%Trsp      = vegn%dailyTrsp
      out_daily_tile%Evap      = vegn%dailyEvap
      out_daily_tile%Runoff    = vegn%dailyRoff
      out_daily_tile%ws1       = vegn%wcl(1)*thksl(1)*1000.
      out_daily_tile%ws2       = vegn%wcl(2)*thksl(2)*1000.
      out_daily_tile%ws3       = vegn%wcl(3)*thksl(3)*1000.
      out_daily_tile%LAI       = vegn%LAI
      out_daily_tile%GPP       = vegn%dailyGPP
      out_daily_tile%Rauto     = vegn%dailyResp
      out_daily_tile%Rh        = vegn%dailyRh
      out_daily_tile%NSC       = vegn%NSC
      out_daily_tile%seedC     = vegn%SeedC
      out_daily_tile%leafC     = vegn%leafC
      out_daily_tile%rootC     = vegn%rootC
      out_daily_tile%SW_C      = vegn%SapwoodC
      out_daily_tile%HW_C      = vegn%woodC
      out_daily_tile%NSN       = vegn%NSN*1000
      out_daily_tile%seedN     = vegn%SeedN*1000
      out_daily_tile%leafN     = vegn%leafN*1000
      out_daily_tile%rootN     = vegn%rootN*1000
      out_daily_tile%SW_N      = vegn%SapwoodN *1000
      out_daily_tile%HW_N      = vegn%WoodN *1000
      out_daily_tile%McrbC     = vegn%MicrobialC
      out_daily_tile%fastSOM   = vegn%metabolicL
      out_daily_tile%slowSOM   = vegn%structuralL
      out_daily_tile%McrbN     = vegn%MicrobialN*1000
      out_daily_tile%fastSoilN = vegn%metabolicN*1000
      out_daily_tile%slowSoilN = vegn%structuralN*1000
      out_daily_tile%mineralN  = vegn%mineralN*1000
      out_daily_tile%N_uptk    = vegn%dailyNup*1000

    endif

    !annual tile
    ! Annual summary:
    vegn%annualNup  = vegn%annualNup  + vegn%dailyNup
    vegn%annualGPP  = vegn%annualGPP  + vegn%dailygpp
    vegn%annualNPP  = vegn%annualNPP  + vegn%dailynpp
    vegn%annualResp = vegn%annualResp + vegn%dailyresp
    vegn%annualRh   = vegn%annualRh   + vegn%dailyrh
    vegn%annualPrcp = vegn%annualPrcp + vegn%dailyPrcp
    vegn%annualTrsp = vegn%annualTrsp + vegn%dailytrsp
    vegn%annualEvap = vegn%annualEvap + vegn%dailyevap
    vegn%annualRoff = vegn%annualRoff + vegn%dailyRoff

    ! zero:
    vegn%dailyNup  = 0.0
    vegn%dailyGPP  = 0.0
    vegn%dailyNPP  = 0.0
    vegn%dailyResp = 0.0
    vegn%dailyRh   = 0.0
    vegn%dailyPrcp = 0.0
    vegn%dailyTrsp = 0.0
    vegn%dailyEvap = 0.0
    vegn%dailyRoff = 0.0

  end subroutine daily_diagnostics

!======================================================
  subroutine annual_diagnostics(vegn, iyears, out_annual_cohorts, out_annual_tile)
   
   use md_interface_lm3ppa, only: outtype_annual_cohorts, outtype_annual_tile, myinterface

   type(vegn_tile_type), intent(inout) :: vegn
   integer, intent(in) :: iyears
   type(outtype_annual_cohorts), dimension(out_max_cohorts) :: out_annual_cohorts
   type(outtype_annual_tile) :: out_annual_tile
   ! type(spec_data_type) :: sp

!   --------local var --------
    type(cohort_type), pointer :: cc
    real treeG, fseed, fleaf, froot,fwood,dDBH
    real :: plantC, plantN, soilC, soilN
    integer :: i

    ! re-initialise to avoid elements not updated when number 
    ! of cohorts declines from one year to the next
    out_annual_cohorts(:)%year    = dummy
    out_annual_cohorts(:)%cID     = dummy
    out_annual_cohorts(:)%PFT     = dummy
    out_annual_cohorts(:)%layer   = dummy
    out_annual_cohorts(:)%density = dummy
    out_annual_cohorts(:)%f_layer = dummy
    out_annual_cohorts(:)%dDBH    = dummy
    out_annual_cohorts(:)%dbh     = dummy
    out_annual_cohorts(:)%height  = dummy
    out_annual_cohorts(:)%age     = dummy
    out_annual_cohorts(:)%Acrown  = dummy
    out_annual_cohorts(:)%wood    = dummy
    out_annual_cohorts(:)%nsc     = dummy
    out_annual_cohorts(:)%NSN     = dummy
    out_annual_cohorts(:)%NPPtr   = dummy
    out_annual_cohorts(:)%seed    = dummy
    out_annual_cohorts(:)%NPPL    = dummy
    out_annual_cohorts(:)%NPPR    = dummy
    out_annual_cohorts(:)%NPPW    = dummy
    out_annual_cohorts(:)%GPP     = dummy
    out_annual_cohorts(:)%NPP     = dummy
    out_annual_cohorts(:)%N_uptk  = dummy
    out_annual_cohorts(:)%N_fix   = dummy
    out_annual_cohorts(:)%maxLAI  = dummy
    out_annual_cohorts(:)%Volume  = dummy

    ! Cohotrs ouput
    do i = 1, vegn%n_cohorts

      cc => vegn%cohorts(i)
      treeG = cc%seedC + cc%NPPleaf + cc%NPProot + cc%NPPwood
      fseed = cc%seedC/treeG
      fleaf = cc%NPPleaf/treeG
      froot = cc%NPProot/treeG
      fwood = cc%NPPwood/treeG
      dDBH  = (cc%DBH - cc%DBH_ys)*1000
      cc%Volume = (cc%bsw+cc%bHW)/spdata(cc%species)%rho_wood

      out_annual_cohorts(i)%year    = iyears
      out_annual_cohorts(i)%cID     = cc%ccID
      out_annual_cohorts(i)%PFT     = cc%species
      out_annual_cohorts(i)%layer   = cc%layer
      out_annual_cohorts(i)%density = cc%nindivs*10000
      out_annual_cohorts(i)%f_layer = cc%layerfrac
      out_annual_cohorts(i)%dDBH    = dDBH
      out_annual_cohorts(i)%dbh     = cc%dbh
      out_annual_cohorts(i)%height  = cc%height
      out_annual_cohorts(i)%age     = cc%age
      out_annual_cohorts(i)%Acrown  = cc%crownarea
      out_annual_cohorts(i)%wood    = cc%bsw+cc%bHW
      out_annual_cohorts(i)%nsc     = cc%nsc
      out_annual_cohorts(i)%NSN     = cc%NSN*1000
      out_annual_cohorts(i)%NPPtr   = treeG
      out_annual_cohorts(i)%seed    = fseed
      out_annual_cohorts(i)%NPPL    = fleaf
      out_annual_cohorts(i)%NPPR    = froot
      out_annual_cohorts(i)%NPPW    = fwood
      out_annual_cohorts(i)%GPP     = cc%annualGPP
      out_annual_cohorts(i)%NPP     = cc%annualNPP
      out_annual_cohorts(i)%N_uptk  = cc%annualNup*1000
      out_annual_cohorts(i)%N_fix   = cc%annualfixedN*1000
      out_annual_cohorts(i)%maxLAI  = spdata(cc%species)%laimax
      out_annual_cohorts(i)%Volume  = cc%Volume

    enddo

    ! tile pools output
    call summarize_tile( vegn )

    do i = 1, vegn%n_cohorts
      cc => vegn%cohorts(i)
      vegn%annualfixedN = vegn%annualfixedN  + cc%annualfixedN * cc%nindivs
      
    enddo

    plantC    = vegn%NSC + vegn%SeedC + vegn%leafC + vegn%rootC + vegn%SapwoodC + vegn%woodC
    plantN    = vegn%NSN + vegn%SeedN + vegn%leafN + vegn%rootN + vegn%SapwoodN + vegn%woodN
    soilC     = vegn%MicrobialC + vegn%metabolicL + vegn%structuralL
    soilN     = vegn%MicrobialN + vegn%metabolicN + vegn%structuralN + vegn%mineralN
    vegn%totN = plantN + soilN

    out_annual_tile%year       = iyears
    out_annual_tile%CAI        = vegn%CAI
    out_annual_tile%LAI        = vegn%LAI
    out_annual_tile%density    = vegn%nindivs*10000 !xxx New tile out
    out_annual_tile%DBH        = vegn%DBH
    out_annual_tile%density12  = vegn%nindivs12*10000 !xxx New tile out
    out_annual_tile%DBH12      = vegn%DBH12
    out_annual_tile%QMD        = vegn%QMD
    out_annual_tile%NPP        = vegn%annualNPP !xxx New tile out
    out_annual_tile%GPP        = vegn%annualGPP
    out_annual_tile%Rauto      = vegn%annualResp
    out_annual_tile%Rh         = vegn%annualRh
    out_annual_tile%rain       = vegn%annualPrcp
    out_annual_tile%SoilWater  = vegn%SoilWater
    out_annual_tile%Transp     = vegn%annualTrsp
    out_annual_tile%Evap       = vegn%annualEvap
    out_annual_tile%Runoff     = vegn%annualRoff
    out_annual_tile%plantC     = plantC
    out_annual_tile%soilC      = soilC
    out_annual_tile%plantN     = plantN *1000
    out_annual_tile%soilN      = soilN * 1000
    out_annual_tile%totN       = (plantN+soilN)*1000
    out_annual_tile%NSC        = vegn%NSC
    out_annual_tile%SeedC      = vegn%SeedC
    out_annual_tile%leafC      = vegn%leafC
    out_annual_tile%rootC      = vegn%rootC
    out_annual_tile%SapwoodC   = vegn%SapwoodC
    out_annual_tile%WoodC      = vegn%woodC
    out_annual_tile%NSN        = vegn%NSN*1000
    out_annual_tile%SeedN      = vegn%SeedN*1000
    out_annual_tile%leafN      = vegn%leafN*1000
    out_annual_tile%rootN      = vegn%rootN*1000
    out_annual_tile%SapwoodN   = vegn%SapwoodN *1000
    out_annual_tile%WoodN      = vegn%WoodN *1000
    out_annual_tile%McrbC      = vegn%MicrobialC
    out_annual_tile%fastSOM    = vegn%metabolicL
    out_annual_tile%SlowSOM    = vegn%structuralL
    out_annual_tile%McrbN      = vegn%MicrobialN*1000
    out_annual_tile%fastSoilN  = vegn%metabolicN*1000
    out_annual_tile%slowSoilN  = vegn%structuralN*1000
    out_annual_tile%mineralN   = vegn%mineralN*1000
    out_annual_tile%N_fxed     = vegn%annualfixedN*1000
    out_annual_tile%N_uptk     = vegn%annualNup*1000
    out_annual_tile%N_yrMin    = vegn%annualN*1000
    out_annual_tile%N_P2S      = vegn%N_P2S_yr*1000
    out_annual_tile%N_loss     = vegn%Nloss_yr*1000
    out_annual_tile%totseedC   = vegn%totseedC*1000
    out_annual_tile%totseedN   = vegn%totseedN*1000
    out_annual_tile%Seedling_C = vegn%totNewCC*1000
    out_annual_tile%Seedling_N = vegn%totNewCN*1000
    out_annual_tile%MaxAge     = vegn%MaxAge
    out_annual_tile%MaxVolume  = vegn%MaxVolume
    out_annual_tile%MaxDBH     = vegn%MaxDBH

    ! I cannot figure out why N losing. Hack!
   if(myinterface%params_siml%do_closedN_run) call Recover_N_balance(vegn)

 end subroutine annual_diagnostics

subroutine Recover_N_balance(vegn)

   type(vegn_tile_type), intent(inout) :: vegn

      if(abs(vegn%totN-vegn%initialN0)*1000>0.001)then
         vegn%structuralN = vegn%structuralN - vegn%totN + vegn%initialN0
         vegn%totN =  vegn%initialN0
      endif

 end subroutine

end module datatypes
