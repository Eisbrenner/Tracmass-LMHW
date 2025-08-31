SUBROUTINE setup_grid
  !* =============================================================
  !* Set up the grid
  !* =============================================================
  !* Subroutine for defining the grid of the GCM. Run once
  !* before the loop starts.
  !* -------------------------------------------------------------
  !* The following arrays have to be populated:
  !*
  !*  dxdy - Area of horizontal cell-walls.
  !*  dzt  - Height of k-cells in 4 dim
  !*  kmt  - Number of k-cells from surface to seafloor.
  !*
  !* The following might be needed to calculate
  !* dxdy, uflux, and vflux
  !*
  !*  dzu - Height of each u-gridcell.
  !*  dzv - Height of each v-gridcell.
  !*  dxv -
  !*  dyu -
  !* -------------------------------------------------------------

  USE mod_precdef
  USE mod_param

  USE mod_grid

  IMPLICIT none

  INTEGER               :: jj, ii, kk
  REAL(PP)              :: dlon, dlat, Tminlat
  REAL(PP), ALLOCATABLE :: dz(:)

  !* Minimum latitude of the T grid
  Tminlat = 9.875d0

  !* Twelvs degree resolution
  dlon = 1./12.; dlat = 1./12.

  !* dx and dy in u and v points & grid area
  dxv(:, :) = 0.d0; dyu(:, :) = 0.d0

  DO jj = 1, jmt
     DO ii = 1, imt
        dx = dlon*deg*COS((Tminlat+0.5*dlat+dlat*(jj - 0.5))*radian)
        dy = dlat*deg

        dxv(ii, jj) = dlon*deg*COS((Tminlat+0.5*dlat+dlat*jj)*radian)
        dyu(ii, jj) = dy

        dxdy(ii, jj) = dx*dy

     END DO
  END DO

  !* Read in dz
  ALLOCATE (dz(1:KM))

  OPEN (12, FILE=trim(topoDataDir)//trim(zgridFile))

  DO kk = 1, KM
     READ (12, *) dz(kk)
  END DO

  CLOSE (12)

  do kk = 1, KM
     dzt(:, :, kk, :) = dz(kk)
  end do

END SUBROUTINE setup_grid
