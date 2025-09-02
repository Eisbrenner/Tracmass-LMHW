SUBROUTINE read_field
   !*==========================================================================
   !*
   !* Purpose
   !* -------
   !*
   !* Read GLORYS12 model output to advect trajectories.
   !* Will be called by the loop each time step.
   !*
   !* Method
   !* ------
   !*
   !* Read velocities and optionally some tracers from netCDF files and
   !* update velocity fields for TRACMASS.
   !*
   !* Updates the variables:
   !*   uflux and vflux
   !* ==========================================================================

   USE mod_precdef
   USE mod_param
   USE mod_vel
   USE mod_time
   USE mod_grid
   USE mod_getfile
   USE mod_tracervars
   USE mod_tracers
   USE mod_calendar
   USE mod_swap

   USE netcdf

   IMPLICIT none

   INTEGER        :: kk, itrac, i0, j0, ios, jos

   REAL(DP), ALLOCATABLE, DIMENSION(:, :, :)  :: tmpu, tmpv
   CHARACTER(len=200)                         :: fieldFile, dateprefix

   !* Reassign the time index of uflux and vflux, dzt, dzdt, hs, ...
   CALL swap_time()

   !* Data files
   dateprefix = ' '

   ALLOCATE (tmpu(0:imt, 1:jmt, 1:km), tmpv(1:imt, 0:jmt, 1:km))

   !* Reading 3-time step variables
   !* ===========================================================================
   IF (ints == 0) THEN

      !* 1 - Past
      IF (loopYears) THEN
         !! removed, only needed for loopYears == True
         STOP "loopYears not implemented."
      END IF

      !* 2 - Present
      ! time position
      nctstep = currDay

      dateprefix = filledFileName(dateFormat, currYear, currMon, currDay)

      !* Tracers
      IF (l_tracers) THEN
         !! ignore this for now
         STOP "l_tracers not implemented."
      END IF

   END IF

   IF (ints < intrun - 1 .OR. loopYears) THEN
      !* 3 - Future
      ! time position
      ! this is 1 because I only have one time step per file
      nctstep = nextDay

      dateprefix = filledFileName(dateFormat, nextYear, nextMon, nextDay)

      !* Tracers
      IF (l_tracers) THEN
         !! ignore this for now
         STOP "l_tracers not implemented."
      END IF

   END IF

   !* Reading 2-time step variables
   !* In this case: velocities and tracers
   !* ===========================================================================

   nctstep = currDay

   dateprefix = filledFileName(dateFormat, currYear, currMon, currDay)

   fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(dateprefix)//TRIM(uGridName)//TRIM(fileSuffix)
   tmpu(0:imt, 1:jmt, 1:km) = get3DfieldNC(fieldFile, ueul_name, [imindom, jmindom, 1, nctstep], [imt + 1, jmt, km, 1], 'st')

   fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(dateprefix)//TRIM(vGridName)//TRIM(fileSuffix)
   tmpv(1:imt, 0:jmt, 1:km) = get3DfieldNC(fieldFile, veul_name, [imindom, jmindom, 1, nctstep], [imt, jmt + 1, km, 1], 'st')

   !* Tracers
   IF (l_tracers) THEN
      !! ignore this for now
      STOP "l_tracers not implemented."
   END IF

   ! uflux and vflux computation
   FORALL (kk=1:km) uflux(1:imt, :, kk, 2) = tmpu(1:imt, :, kk)*dyu(:, :)*dzt(:, :, kk, 2)
   FORALL (kk=1:km) vflux(:, 1:jmt, kk, 2) = tmpv(:, 1:jmt, kk)*dxv(:, :)*dzt(:, :, kk, 2)

   !* dzdt calculation
   dzdt = 0.d0

   !* Reverse the sign of fluxes if trajectories are run backward in time.
   CALL swap_sign()

END SUBROUTINE read_field
