! example of logger subroutine in salmon !uemoto!
subroutine salmon_open_logfile() !uemoto!
use salmon_parallel, only: nproc_id_global !uemoto!
implicit none !uemoto!
character(100) :: logfile !uemoto!
write(logfile,'(i3.3, ".log")') nproc_id_global !uemoto!
open(unit=777, file=trim(logfile)) !uemoto!
return !uemoto!
end subroutine salmon_open_logfile !uemoto!
program main
  use salmon_global
  use salmon_parallel
  use inputoutput
  use math_constants
  implicit none

  call set_math_constants

  call setup_parallel
  if (nproc_id_global == 0) then
    call print_software_version
  endif

  call read_input
  call  salmon_open_logfile() !uemoto!

  select case(iperiodic)
  case(0)
    call gceed
  case(3)
    select case(domain_parallel)
    case('y')
      call gceed
    case('n')
      call arted
    case default
      stop 'invalid domain_parallel'
    end select
  case default
    stop 'invalid iperiodic'
  end select

  call end_parallel
contains
  subroutine print_software_version
    use salmon_xc, only: print_xc_info
    implicit none
    include 'versionf.h'
    print '(A)',         '##############################################################################'
    print '(A)',         '# SALMON: Scalable Ab-initio Light-Matter simulator for Optics and Nanoscience'
    print '(A)',         '#'
    print '(A,I1,".",I1,".",I1)', &
    &                    '#                             Version ', SALMON_VER_MAJOR, SALMON_VER_MINOR, SALMON_VER_MICRO
    if (GIT_FOUND) then 
      print '(A)',       '#'
      print '(A,A,A,A)', '#   [Git revision] ', GIT_COMMIT_HASH, ' in ', GIT_BRANCH
    endif
    print '(A)',         '##############################################################################'
    
    call print_xc_info()    
  end subroutine
end program main
