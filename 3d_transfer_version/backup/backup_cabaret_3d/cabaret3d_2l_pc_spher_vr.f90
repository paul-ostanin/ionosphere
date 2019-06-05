! Solve mass continuity equation on the sphere
! 3d spherical layer case, C-grid,
! uniform grid in longitude,latitude and nonuniform in height

program cabaret2d

use transfer_mod, only : init_transfer, step_of_transfer

implicit none

call init_transfer()

call step_of_transfer()


end
