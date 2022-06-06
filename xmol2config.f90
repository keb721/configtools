!---------------------------------------------------------------------------!
!  Reads specified xmol file and prints equivalent DL_POLY CONFIG file to   !
!  to stdout. Expects to find the matrix of cell vectors on the comment     !
!  line of the input xmol file, listed as 9 components.                     !
!---------------------------------------------------------------------------!
! version 1.0 - D.Quigley ( University of Warwick )                         !
!---------------------------------------------------------------------------!
program xmol2config

  implicit none 
  integer,parameter :: dp = kind(1.0d0)

  integer :: N,idim,ierr,iatom,num_args,iarg

  real(kind=dp),dimension(3,3) :: Hmatrix
  real(kind=dp) :: xxx,yyy,zzz

  character(30), dimension(0:10) :: command_line
  character(8)   :: dumchar
  character(30)  :: filename

  integer,parameter :: xml = 25 ! unit number for input xmol file
  integer,parameter :: cfg = 26 ! unit number for output CONFIG
  

  ! check that there is one arguments.
  num_args = iargc()

  if (num_args/=1) then
     write(0,*)
     write(0,*) '          X M O L 2 C O N F I G             '
     write(0,*)
     write(0,*) '    Usage: xmol2config <xmolfile>           '
     write(0,*)
     write(0,*) '    D. Quigley - University of Warwick      '
     write(0,*)
     stop
  end if

  ! get command line arguments
  do iarg = 1, num_args
     call getarg(iarg,command_line(iarg))
  end do

  filename = trim(command_line(1))
  
  ! open the specified xmol file
  open(unit=xml,file=trim(filename),status='old',iostat=ierr)
  if (ierr/=0) stop 'Error opening specified xmol file for input'

  ! open the output file
  open(unit=cfg,file='CONFIG',status='replace',iostat=ierr)
  if (ierr/=0) stop 'Error opening CONFIG for output'

  ! read the number of atoms
  read(xml,*,iostat=ierr)N
  if (ierr/=0) stop 'Error reading number of atoms on line 1 of config.xmol'

  ! read the matrix of cell vectors
  ! sometimes there's a character before the cell matrix
  read(xml,*,iostat=ierr) dumchar, hmatrix           ! blank comment line in xmol format
  if (ierr/=0) then
     rewind(xml)
     read(xml,*)
     read(xml,*,iostat=ierr) hmatrix
     if (ierr/=0) then
        stop 'Could not extract matrix of cell vectors from line 2 of xmolfile'
     end if
  end if
  
  ! write the header of the config file
  write(cfg,*)"CONFIG file created from Xmol file "//trim(filename)
  write(cfg,'(2i10)')0,3
  
  ! write the cell vector information to the CONFIG file
  do idim = 1,3
     write(cfg,'(3f20.6)')hmatrix(:,idim)
  end do

  ! loop over the atoms in the file
  do iatom = 1,N

     read(xml,*)dumchar,xxx,yyy,zzz

     write(cfg,'(a8,i10)')dumchar,iatom
     write(cfg,'(3f20.6)')xxx,yyy,zzz

  end do

  close(xml)
  close(cfg)

end program xmol2config
