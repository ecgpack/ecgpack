module wp_def
!This module contains basic constants and basic input/output subroutines
!which depend on the choice of precision, compiler, and MPI implementation

!include 'mpif.h'
  use mpi

!Define kind parameter for real type
  integer,parameter     :: wp=8

!This is data type identifier for MPI corresponding to real type of kind wp
  integer,parameter    :: MPI_WP=MPI_DOUBLE_PRECISION

!This is the number of particles in the system that should be set by the user.
!For reasons related to performance, it is made a fixed (compile time) parameter.
  integer,parameter :: Glob_AllowedNumOfParticles=4

contains

!Subroutines writereal and writerealarr (writerealadv and writerealarradv) realize nonadvanced
!(advanced) output of real(wp) type. For portability and easiness of modification of the
!program, all output of real(wp) type, both on screen and to external files should be done
!via calling these subroutines.

  subroutine writereal(u,r)
    integer u          !i/o unit
    real(wp) r      !real number that needs to be written
    write(u,'(1x,e23.16)',advance='no') r
  end subroutine writereal

  subroutine writerealadv(u,r)
    integer u          !i/o unit
    real(wp) r      !real number that needs to be written
    write(u,'(1x,e23.16)') r
  end subroutine writerealadv

  subroutine writerealarr(u,r,k)
    integer u          !i/o unit
    real(wp) r(k)   !real array that needs to be written
    integer k          !the number of elements to write (writing begins with element 1)
    integer i
    do i=1,k
      write(u,'(1x,e23.16)',advance='no') r(i)
    enddo
  end subroutine writerealarr

  subroutine writerealarradv(u,r,k)
    integer u          !i/o unit
    real(wp) r(k)   !real array that needs to be written
    integer k          !the number of elements to write (writing begins with element 1)
    integer i
    do i=1,k-1
      write(u,'(1x,e23.16)',advance='no') r(i)
    enddo
    write(u,'(1x,e23.16)') r(k)
  end subroutine writerealarradv

!Subroutines writestring and writestringadv realize nonadvanced (advanced)
!output of strings

  subroutine writestring(u,s,k)
    integer u          !i/o unit
    character(*) s     !string that needs to be written
    integer k          !the length of the string (k first characters) to write.
    integer i
    write(u,'(1x)',advance='no')
    do i=1,k
      write(u,'(a1)',advance='no') s(i:i)
    enddo
  end subroutine writestring

  subroutine writestringadv(u,s,k)
    integer u          !i/o unit
    character(*) s     !string that needs to be written
    integer k          !the length of the string (k first characters) to write.
    integer i
    write(u,'(1x)',advance='no')
    do i=1,k-1
      write(u,'(a1)',advance='no') s(i:i)
    enddo
    write(u,'(a1)') s(k:k)
  end subroutine writestringadv

end module wp_def
