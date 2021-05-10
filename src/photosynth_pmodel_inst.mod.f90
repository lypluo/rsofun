module md_photosynth_inst

  use md_photosynth, only: pmodel, zero_pmodel, outtype_pmodel, calc_ftemp_inst_vcmax, calc_ftemp_inst_jmax, &
                           co2_to_ca, calc_gammastar, calc_kmm
  
  implicit none 
  
  private
  public pmodel_inst, outtype_pmodel_inst
  
  type outtype_pmodel_inst
      ! Returned variables
      real :: assim
      real :: vcmax
      real :: jmax
      ! real :: rd
      ! real :: ci
      ! real :: a_j
      ! real :: a_c
      ! real :: ci_j
      ! real :: ci_c
      ! real :: gpp
  end type outtype_pmodel_inst
  
  
  contains
    
    function pmodel_inst(vcmax25, jmax25, xi, tc, vpd, co2, fapar, ppfd, patm, kphio, tc_home, fixedCi) result(out_pmodel_inst)
      
      !--------------------------------!
      ! Definition of variables        !
      !--------------------------------!
  
      ! Input variables
      real, intent(in) :: vcmax25
      real, intent(in) :: jmax25
      real, intent(in) :: xi
      real, intent(in) :: tc
      real, intent(in) :: vpd
      real, intent(in) :: co2
      real, intent(in) :: fapar
      real, intent(in) :: ppfd
      real, intent(in) :: patm
      real, intent(in) :: kphio
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
  
      ! Output variables
      type(outtype_pmodel_inst) :: out_pmodel_inst
  
      !--------------------------------!
      ! Calcuation of instant response !
      !--------------------------------!
  
      ! Ambient CO2
      ca = co2_to_ca(co2, patm)
      
      ! Instantaneous gammastar
      gammastar = calc_gammastar(tc, patm)
      
      ! Instantaneous K
      kmm = calc_kmm(tc, patm)
        
      ! Instantaneous vcmax and jmax
      vcmax = calc_ftemp_inst_vcmax(tc, tc, tcref = 25.0 ) * vcmax25
      jmax  = calc_ftemp_inst_jmax(tc, tc, tc_home, tcref = 25.0 ) * jmax25
        
      ! Instantaneous kphio - TODO: Not needed because kphio does not adapt on daily changes?
      ! kphio = kphio * calc_ftemp_kphio()
        
      ! Aj, gs free
      L = 1.0 / sqrt(1.0 + ((4.0 * kphio * ppfd) / jmax)**2)
      kv = (ca - gammastar) / (1 + xi / sqrt(vpd))
      ci_j = ca - kv

      ! If check for fixedCi
      if (fixedCi) then
          ci_j = 275e-06 !275ppm
      end if

      a_j = L * kphio * ppfd * (ci_j - gammastar)/(ci_j + 2 * gammastar)
      gs_j = a_j / kv
      
      ! If check for fixedCi
      ! Ac, gs free
      ci_c = ci_j

      if (fixedCi) then
          ci_c = 275e-06 !275ppm
      end if

      a_c = vcmax * (ci_c - gammastar)/(ci_c + kmm)
      gs_c = a_j / kv
        
      ! TODO: rpmodel has unfinished tryout-code here for gs being not free
        
      ! Assimilation
      assim = min(a_j, a_c)
      ci = max(ci_c, ci_j)
  

      !--------------------------------!
      ! Definition of output           !
      !--------------------------------!
      out_pmodel_inst%assim = assim
      out_pmodel_inst%vcmax = vcmax
      out_pmodel_inst%jmax  = jmax
  
    end function pmodel_inst
end module md_photosynth_inst
