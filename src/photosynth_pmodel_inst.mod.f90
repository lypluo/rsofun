module md_photosynth_inst

  use md_photosynth, only: pmodel, zero_pmodel, outtype_pmodel, calc_ftemp_inst_vcmax, calc_ftemp_inst_jmax, co2_to_ca, calc_gammastar, calc_kmm
  
  implicit none 
  
  private
  public pmodel_inst, outtype_pmodel_inst
  
  type outtype_pmodel_inst
      ! Returned variables
      real :: assim
      real :: vcmax
      real :: jmax
      real :: rd
      real :: ci
      real :: gammastar
      real :: anet
      ! real :: a_j
      ! real :: a_c
      ! real :: ci_j
      ! real :: ci_c
      ! real :: gpp

      ! TODO: For debugging
      real :: term

  end type outtype_pmodel_inst
  
  
  contains
    
    function pmodel_inst(vcmax25, jmax25, xi, tc_leaf, vpd, co2, fapar, ppfd, patm, kphio, tc_growth, tc_home, fixedCi) result(out_pmodel_inst)
      
      !--------------------------------!
      ! Definition of variables        !
      !--------------------------------!
  
      ! Input variables
      real, intent(in) :: vcmax25
      real, intent(in) :: jmax25
      real, intent(in) :: xi
      real, intent(in) :: tc_leaf
      real, intent(in) :: vpd
      real, intent(in) :: co2
      real, intent(in) :: fapar
      real, intent(in) :: ppfd
      real, intent(in) :: patm
      real, intent(in) :: kphio
      real, intent(in) :: tc_growth
      real, intent(in) :: tc_home
      logical, intent(in) :: fixedCi
  
      ! Local variables
      real :: ca
      real :: gammastar
      real :: kmm
      real :: jmax
      real :: vcmax
      real :: L
      real :: kv
      real :: ci_j
      real :: ci_c
      real :: ci
      real :: a_j
      real :: a_c
      real :: gs_j
      real :: gs_c
      real :: assim
      real :: anet
      real :: rd_to_vcmax
      real :: rd
      real :: rd_25
      real :: q10

      ! TODO: For debugging
      real :: term
  
      ! Output variables
      type(outtype_pmodel_inst) :: out_pmodel_inst
  
      !--------------------------------!
      ! Calcuation of instant response !
      !--------------------------------!
  
      ! Ambient CO2
      ca = co2_to_ca(co2, patm)
      
      ! Instantaneous gammastar
      gammastar = calc_gammastar(tc_leaf, patm)
      
      ! Instantaneous K
      kmm = calc_kmm(tc_leaf, patm)
        
      ! Instantaneous vcmax and jmax
      ! calc_ftemp_inst_vcmax( tcleaf, tcgrowth, tcref )
      ! calc_ftemp_inst_jmax( tcleaf, tcgrowth, tc_home, tcref )

      vcmax = calc_ftemp_inst_vcmax(tc_leaf, tc_growth, tcref = 25.0 )         * vcmax25
      jmax  = calc_ftemp_inst_jmax(tc_leaf, tc_growth, tc_home, tcref = 25.0 ) * jmax25
        
      ! Aj, gs free
      L = 1.0 / sqrt(1.0 + ((4.0 * kphio * ppfd) / jmax)**2)
      kv = (ca - gammastar) / (1 + xi / sqrt(vpd))
      ci_j = ca - kv

      if (fixedCi) then ! If check for fixedCi
          ci_j = 27.5   ! Corresponds to 275 ppm compared to ambient 10^6 Pa
      end if

      a_j = L * kphio * ppfd * (ci_j / (ci_j + 2*gammastar)) * (1.0 - gammastar/ci_j)


      if (.false.) then ! Taking Jmax limitation from Farquhar (1989) with a curvature parameter of 0.85
        a_j = (kphio * ppfd + jmax - sqrt(( kphio * ppfd + jmax)**2 - (4*kphio*0.85*ppfd*jmax))) / (2*0.85) / 4 * (ci_c - gammastar)/(ci_c + kmm) * (1.0 - gammastar/ci_c) 
      end if

      gs_j = a_j / kv
      
      ! Ac, gs free
      ci_c = ci_j
      a_c = vcmax * (ci_c - gammastar)/(ci_c + kmm)
      gs_c = a_j / kv
                
      ! Assimilation
      assim = min(a_j, a_c)
      ci = max(ci_c, ci_j)

      ! Dark respiration
      rd_to_vcmax = 0.015                   ! Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous TODO: rd_to_vcmax clean-up
      rd = vcmax25 * rd_to_vcmax * calc_ftemp_inst_rd(tc_leaf)

      ! Net Assimilation
      anet = assim - rd

      !--------------------------------!
      ! Definition of output           !
      !--------------------------------!
      out_pmodel_inst%anet  = anet
      out_pmodel_inst%rd    = rd
      out_pmodel_inst%vcmax = vcmax
      out_pmodel_inst%jmax  = jmax
      out_pmodel_inst%gammastar  = gammastar
      out_pmodel_inst%ci    = ci

      ! TODO: For debugging
      out_pmodel_inst%assim = assim
      out_pmodel_inst%rd    = rd
      out_pmodel_inst%term  = (1.0 - gammastar/ci)

    end function pmodel_inst

    ! TODO: Clean up where calc_ftemp_inst_rd() is located and called!

    function calc_ftemp_inst_rd( tc ) result( fr )
    !-----------------------------------------------------------------------
    ! Output:   Factor fr to correct for instantaneous temperature response
    !           of Rd (dark respiration) for:
    !
    !               Rd(temp) = fr * Rd(25 deg C) 
    !
    ! Ref:      Heskel et al. (2016) used by Wang Han et al. (in prep.)
    !-----------------------------------------------------------------------
    ! arguments
    real, intent(in) :: tc      ! temperature (degrees C)

    ! function return variable
    real :: fr                  ! temperature response factor, relative to 25 deg C.

    ! loal parameters
    real, parameter :: apar = 0.1012
    real, parameter :: bpar = 0.0005
    real, parameter :: tk25 = 298.15 ! 25 deg C in Kelvin

    ! local variables
    real :: tk                  ! temperature (Kelvin)

    ! conversion of temperature to Kelvin
    tk = tc + 273.15

    fr = exp( apar * (tc - 25.0) - bpar * (tc**2 - 25.0**2) )
    
  end function calc_ftemp_inst_rd
end module md_photosynth_inst
