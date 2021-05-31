module md_biosphere_lm3ppa

  use datatypes
  use md_vegetation_lm3ppa
  use md_soil_lm3ppa
  use md_params_core

  implicit none

  private
  public biosphere_annual

   type(vegn_tile_type),  pointer :: vegn   
   type(soil_tile_type),  pointer :: soil
   type(cohort_type),     pointer :: cx, cc

contains

  subroutine biosphere_annual(out_biosphere)
    !////////////////////////////////////////////////////////////////
    ! function BIOSPHERE_annual calculates net ecosystem exchange (nee)
    ! in response to environmental boundary conditions (atmospheric 
    ! CO2, temperature, Nitrogen deposition. This SR "replaces" 
    ! LPJ, also formulated as subroutine.
    ! Copyright (C) 2015, see LICENSE, Benjamin David Stocker
    ! contact: b.stocker@imperial.ac.uk
    !----------------------------------------------------------------
    use md_interface_lm3ppa, only: myinterface, outtype_biosphere  
    ! return variable
    type(outtype_biosphere) :: out_biosphere

    ! ! local variables
    integer :: dm, moy, doy
    logical, save :: init = .true.   ! is true only on the first day of the simulation 
    logical, parameter :: verbose = .false.       ! change by hand for debugging etc.

    !----------------------------------------------------------------
    ! Biome-E stuff
    !----------------------------------------------------------------
    integer, parameter :: rand_seed = 86456
    integer, parameter :: totalyears = 10
    integer, parameter :: nCohorts = 1
    real    :: tsoil, soil_theta
    integer :: year0
    integer :: i, j
    integer :: idata
    integer, save :: simu_steps !, datalines
    integer, save :: iyears
    integer, save :: idays
    integer, save :: idoy
    ! character(len=50) :: namelistfile = '~/sofun/params/parameters_Allocation.nml' !'parameters_WC_biodiversity.nml' ! 'parameters_CN.nml'

    !----------------------------------------------------------------
    ! INITIALISATIONS
    !----------------------------------------------------------------
    if (myinterface%steering%init) then

      ! Parameter initialization: Initialize PFT parameters
      call initialize_PFT_data()

      ! Initialize vegetation tile and plant cohorts
      allocate( vegn )
      call initialize_vegn_tile( vegn, nCohorts)
      
      ! Sort and relayer cohorts
      call relayer_cohorts( vegn )

      ! initialise outputs 
      call Zero_diagnostics( vegn )

      year0 = myinterface%climate(1)%year  ! forcingData(1)%year
      iyears = 1
      idoy = 0
      idays  = 0

    endif 

    simu_steps = 0

    !----------------------------------------------------------------
    ! LOOP THROUGH MONTHS
    !----------------------------------------------------------------
    doy = 0
    monthloop: do moy=1,nmonth

      !----------------------------------------------------------------
      ! LOOP THROUGH DAYS
      !----------------------------------------------------------------
      dayloop: do dm=1,ndaymonth(moy)
        
        doy = doy + 1
        idoy = idoy + 1

        if (verbose) print*,'----------------------'
        if (verbose) print*,'YEAR, DOY ', myinterface%steering%year, doy
        if (verbose) print*,'----------------------'

        !----------------------------------------------------------------
        ! FAST TIME STEP
        !----------------------------------------------------------------
        ! get daily mean temperature from hourly/half-hourly data
        vegn%Tc_daily = 0.0
        tsoil         = 0.0
        fastloop: do i=1,myinterface%steps_per_day

          idata = simu_steps + 1
          year0 =  myinterface%climate(idata)%year  ! Current year
          vegn%Tc_daily = vegn%Tc_daily + myinterface%climate(idata)%Tair

          tsoil         = myinterface%climate(idata)%tsoil
          simu_steps    = simu_steps + 1

          !----------------------------------------------------------------
          ! Sub-daily time step at resolution given by forcing (can be 1 = daily)
          !----------------------------------------------------------------
          call vegn_CNW_budget( vegn, myinterface%climate(idata), init )
          call hourly_diagnostics( vegn, myinterface%climate(idata), iyears, idoy, i, out_biosphere%hourly_tile(idata) )
          init = .false.

        enddo fastloop ! hourly or half-hourly

        !-------------------------------------------------
        ! Daily calls after fast loop
        !-------------------------------------------------
        vegn%Tc_daily = vegn%Tc_daily/myinterface%steps_per_day
        tsoil         = tsoil/myinterface%steps_per_day
        soil_theta    = vegn%thetaS

        ! sum over fast time steps and cohorts
        call daily_diagnostics( vegn, iyears, idoy, out_biosphere%daily_cohorts(doy,:), out_biosphere%daily_tile(doy) )

        ! Determine start and end of season and maximum leaf (root) mass
        call vegn_phenology( vegn )

        ! Produce new biomass from 'carbon_gain' (is zero afterwards)
        call vegn_growth_EW( vegn )

      end do dayloop

    end do monthloop

    !----------------------------------------------------------------
    ! Annual calls
    !----------------------------------------------------------------
    idoy = 0

    print*,'sim. year  ', iyears
    print*,'real year: ', year0

    if ( myinterface%params_siml%update_annualLAImax ) call vegn_annualLAImax_update( vegn )

    call annual_diagnostics( vegn, iyears, out_biosphere%annual_cohorts(:), out_biosphere%annual_tile)

    !---------------------------------------------
    ! Reproduction and mortality
    !---------------------------------------------        
    ! Kill all individuals in a cohort if NSC falls below critical point
    call vegn_annual_starvation( vegn )

    ! Natural mortality (reducing number of individuals 'nindivs')
    ! (~Eq. 2 in Weng et al., 2015 BG)
    call vegn_nat_mortality( vegn, real( seconds_per_year ))

    ! seed C and germination probability (~Eq. 1 in Weng et al., 2015 BG)
    call vegn_reproduction( vegn )

    !---------------------------------------------
    ! Re-organize cohorts
    !---------------------------------------------
    call kill_lowdensity_cohorts( vegn )
    call relayer_cohorts( vegn )
    call vegn_mergecohorts( vegn )

    ! !---------------------------------------------
    ! ! Set annual variables zero
    ! !---------------------------------------------
    call Zero_diagnostics( vegn )

    ! update the years of model run
    iyears = iyears + 1

    if (myinterface%steering%finalize) then
      !----------------------------------------------------------------
      ! Finazlize run: deallocating memory
      !----------------------------------------------------------------
      deallocate(vegn)

    end if

    ! xxx write restart file all info now in cohort_conserved_type, vegn_tile_conserved_type


    if (verbose) print*,'Done with biosphere for this year. Guete Rutsch!'

  end subroutine biosphere_annual

end module md_biosphere_lm3ppa
