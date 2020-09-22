! This file is used on Cheops and Juwels
! to supplement missing Git/Cmake functionalities
! Variables in this file are filled in by CMake at build time,
! used to record the git version and hash.

module modversion

  implicit none

  character(80) :: git_version='2.17' ! "@GIT_VERSION@"

end module modversion
