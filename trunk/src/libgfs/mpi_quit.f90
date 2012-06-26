subroutine mpi_quit(ierr)
  integer, intent(in) :: ierr
  print *,'stopping with error code,ierr'
  stop 
end subroutine mpi_quit
