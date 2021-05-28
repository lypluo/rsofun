module md_biosphere_pmodel

  use md_params_core
  use md_classdefs
  use md_waterbal, only: waterbal, solar, getpar_modl_waterbal
  use md_gpp_pmodel, only: getpar_modl_gpp, gpp
  use md_vegdynamics_pmodel, only: vegdynamics
  use md_tile_pmodel, only: tile_type, tile_fluxes_type, initglobal_tile, initdaily_tile_fluxes, &
    getpar_modl_tile, diag_daily, diag_annual, init_annual
  use md_plant_pmodel, only: getpar_modl_plant
  ! use md_params_soil_pmodel, only: getpar_soil
  ! use md_soiltemp, only: soiltemp
  use md_sofunutils, only: calc_patm


  implicit none

  private
  public biosphere_annual

  !----------------------------------------------------------------
  ! Module-specific (private) variables
  !----------------------------------------------------------------
  ! derived types from L1 modules
  type(tile_type),        dimension(nlu) :: tile             ! has gridcell-dimension because values are stored between years
  type(tile_fluxes_type), dimension(nlu) :: tile_fluxes      ! has no gridcell-dimension values need not be recorded

contains

  function biosphere_annual() result( out_biosphere )
    !////////////////////////////////////////////////////////////////
    ! function BIOSPHERE_annual calculates net ecosystem exchange (nee)
    ! in response to environmental boundary conditions (atmospheric 
    ! CO2, temperature, Nitrogen deposition. This SR "replaces" 
    ! LPJ, also formulated as subroutine.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_interface_pmodel, only: myinterface, outtype_biosphere
    use md_sofunutils, only: daily2monthly
  
    ! return variable
    type(outtype_biosphere) :: out_biosphere

    ! local variables
    integer :: dm, moy, doy
    logical, save           :: init_daily = .true.   ! is true only on the first day of the simulation 
    logical, parameter      :: verbose = .false.     ! change by hand for debugging etc.

    !----------------------------------------------------------------
    ! INITIALISATIONS
    !----------------------------------------------------------------
    if (myinterface%steering%init) then

      !----------------------------------------------------------------
      ! GET MODEL PARAMETERS
      ! read model parameters that may be varied for optimisation
      !----------------------------------------------------------------
      if (verbose) print*, 'getpar_modl() ...'
      call getpar_modl_tile()
      call getpar_modl_plant()
      call getpar_modl_waterbal()
      call getpar_modl_gpp()
      if (verbose) print*, '... done'

      !----------------------------------------------------------------
      ! Initialise pool variables and/or read from restart file (not implemented)
      !----------------------------------------------------------------
      if (verbose) print*, 'initglobal_() ...'
      call initglobal_tile( tile(:) )
      if (verbose) print*, '... done'

    endif 

    !----------------------------------------------------------------
    ! Set annual sums to zero
    !----------------------------------------------------------------
    call init_annual( tile_fluxes(:) )

    !----------------------------------------------------------------
    ! LOOP THROUGH MONTHS
    !----------------------------------------------------------------
    doy=0
    monthloop: do moy=1,nmonth

      !----------------------------------------------------------------
      ! LOOP THROUGH DAYS
      !----------------------------------------------------------------
      dayloop: do dm=1,ndaymonth(moy)
        doy=doy+1

        if (verbose) print*,'----------------------'
        if (verbose) print*,'YEAR, Doy ', myinterface%steering%year, doy
        if (verbose) print*,'----------------------'

        !----------------------------------------------------------------
        ! initialise daily updated variables 
        !----------------------------------------------------------------
        if (verbose) print*,'calling initdaily_() ...'
        call initdaily_tile_fluxes( tile_fluxes(:) )
        if (verbose) print*,'... done.'

        !----------------------------------------------------------------
        ! Get radiation based on daily temperature, sunshine fraction, and 
        ! elevation.
        !----------------------------------------------------------------
        if (verbose) print*,'calling solar() ... '
        if (verbose) print*,'    with argument lat = ', myinterface%grid%lat
        if (verbose) print*,'    with argument elv = ', myinterface%grid%elv
        if (verbose) print*,'    with argument dfsun (ann. mean) = ', sum( myinterface%climate(:)%dfsun / ndayyear )
        if (verbose) print*,'    with argument dppfd (ann. mean) = ', sum( myinterface%climate(:)%dppfd / ndayyear )
        call solar( tile_fluxes(:), &
                    myinterface%grid, & 
                    myinterface%climate(doy),  &
                    doy &
                    )
        if (verbose) print*,'... done'

        !----------------------------------------------------------------
        ! update canopy and tile variables and simulate daily 
        ! establishment / sprouting
        !----------------------------------------------------------------
        if (verbose) print*,'calling vegdynamics() ... '
        call vegdynamics( tile(:), &
                          myinterface%vegcover(doy)%dfapar, &
                          myinterface%fpc_grid(:) &
                          )
        if (verbose) print*,'... done'

        !----------------------------------------------------------------
        ! calculate GPP
        !----------------------------------------------------------------
        if (verbose) print*,'calling gpp() ... '
        call gpp( tile(:), &
                  tile_fluxes(:), &
                  myinterface%pco2, &
                  myinterface%climate(doy), &
                  myinterface%vegcover(doy), &
                  myinterface%grid, &
                  myinterface%params_siml%soilmstress, &
                  myinterface%params_siml%tempstress, &
                  myinterface%tc_home, &
                  init_daily &
                  )

        if (verbose) print*,'... done'

        !----------------------------------------------------------------
        ! get soil moisture, and runoff
        !----------------------------------------------------------------
        if (verbose) print*,'calling waterbal() ... '
        call waterbal(  tile(:), &
                        tile_fluxes(:), &
                        myinterface%grid, &
                        myinterface%climate(doy), &
                        doy &
                        )
        if (verbose) print*,'... done'

        ! !----------------------------------------------------------------
        ! ! calculate soil temperature
        ! !----------------------------------------------------------------
        ! if (verbose) print*, 'calling soiltemp() ... '
        ! call soiltemp(&
        !               tile(:)%soil, &
        !               myinterface%climate%dtemp(:), &
        !               size(myinterface%grid), &
        !               myinterface%steering%init, &
        !               jpngr, & 
        !               moy, & 
        !               doy & 
        !               )
        ! if (verbose) print*, '... done'

        !----------------------------------------------------------------
        ! daily diagnostics (e.g., sum over plant within canopy)
        !----------------------------------------------------------------
        call diag_daily(tile(:), tile_fluxes(:))

        !----------------------------------------------------------------
        ! populate function return variable
        !----------------------------------------------------------------
        if (nlu>1) stop 'think about nlu > 1'
        
        out_biosphere%fapar(doy)   = tile(1)%canopy%fapar
        out_biosphere%gpp(doy)     = tile_fluxes(1)%canopy%dgpp
        out_biosphere%transp(doy)  = tile_fluxes(1)%canopy%daet
        out_biosphere%latenth(doy) = tile_fluxes(1)%canopy%daet_e
        out_biosphere%pet(doy)     = tile_fluxes(1)%canopy%dpet
        out_biosphere%vcmax(doy)   = tile_fluxes(1)%canopy%vcmax
        out_biosphere%jmax(doy)    = tile_fluxes(1)%canopy%jmax
        out_biosphere%vcmax25(doy) = tile_fluxes(1)%canopy%vcmax25
        out_biosphere%jmax25(doy)  = tile_fluxes(1)%canopy%jmax25
        out_biosphere%gs_accl(doy) = tile_fluxes(1)%canopy%gs_accl
        out_biosphere%wscal(doy)   = tile(1)%soil%phy%wscal

        ! Flexible debugging variables
        out_biosphere%debug1(doy)  = tile_fluxes(1)%canopy%debug1
        out_biosphere%debug2(doy)  = tile_fluxes(1)%canopy%debug2
        out_biosphere%debug3(doy)  = tile_fluxes(1)%canopy%debug3
        out_biosphere%debug4(doy)  = tile_fluxes(1)%canopy%debug4
        out_biosphere%debug5(doy)  = tile_fluxes(1)%canopy%debug5
        out_biosphere%debug6(doy)  = tile_fluxes(1)%canopy%debug6

        init_daily = .false.

      end do dayloop

    end do monthloop

    !----------------------------------------------------------------
    ! annual diagnostics
    !----------------------------------------------------------------
    call diag_annual( tile(:), tile_fluxes(:) )
    

    if (verbose) print*,'Done with biosphere for this year. Guete Rutsch!'

  end function biosphere_annual

end module md_biosphere_pmodel
