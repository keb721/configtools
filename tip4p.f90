!------------------------------------------------------------------------------------!
!   Code to rationalise a water configuration for use with a TIP4P style water       !
!   model. Reads a water configuration from config.xmol and....                      !
!                                                                                    !
!   1. Adjusts O-H bond lengths in that configuration to match that specifed         !
!      in the parameter OH_length.                                                   !
!   2. Adjusts the H-O-H bond angle to match that specified in the parameter         !
!      H-O-H angle.                                                                  !
!   3. Adds the massless charged 'M' site at the specified distance along the        !
!      bisector.                                                                     !
!   4. Unwraps co-ordinates such that no O-H bonds cross a periodic boundary -       !
!      prevents DL_POLY from throwing a major wobbly during quaternion setup.        !
!------------------------------------------------------------------------------------!
! version 1.0 - D.Quigley ( University of Warwick )                                  !  
!------------------------------------------------------------------------------------!
program tip4p

  implicit none

  integer,parameter :: dp = kind(1.0d0)      ! use double precision

  ! Pi 
  real(kind=dp),Parameter :: Pi = 3.141592653589793238462643383279502884197_dp

  ! integers/counters
  integer :: i,l,foundcount,found_this_pass,k,idim,N,j

  ! vectors required for doing the geometry
  real(kind=dp),dimension(3)   :: vec1,vec2,vec3,vec4,cell,oh1,oh2
  real(kind=dp),dimension(3)   :: unit1,unit2

  ! temporary matrix of cell vectors
  real(kind=dp),dimension(3,3) :: dumh

  ! Distance OH common to all TIP4P variants 
  real(kind=dp),parameter :: OH_length = 0.9572_dp

  ! H - O - H bond angle required in radians
  real(kind=dp),parameter :: HOH_angle = 104.52_dp*Pi/180.0_dp
  
  ! Distance OM to use for the current TIP4P variant
  real(kind=dp) :: OM_length = 0.15_dp 

  character(1) :: elemid   ! chemical symbol for the current atom
  character(1) :: dumchar  ! dummy character for cell comment line

  logical :: secondbond    ! is this the second O-H bond for the current H

  real(kind=dp) :: length          ! current bond length
  real(kind=dp) :: theta,costheta  ! current bond angle and cosine thereof
  real(kind=dp) :: mag_oh1,mag_oh2 ! length of O-H bonds before correction

  integer                        :: iarg,idata,ierr=0
  integer                        :: num_args
  character(30), dimension(0:10) :: command_line
  character(30)                  :: filename, style


  ! check that there are two arguments.
  num_args = iargc()

  if (num_args/=2) then
     write(0,*)
     write(0,*) '                T I P 4 P                     '
     write(0,*)
     write(0,*) '    Usage: tip4p <xmolfile> <tip4p_style>     '
     write(0,*)
     write(0,*) '    D. Quigley - University of Warwick        '
     write(0,*)
     stop
  end if

  ! get command line arguments
  do iarg = 1, num_args
     call getarg(iarg,command_line(iarg))
  end do

  filename = trim(command_line(1))
  style    = trim(command_line(2))
  
  ! open specified input file, bail if it's not there
  open(unit=25,file=trim(filename),status='old',iostat=ierr)
  if (ierr/=0) stop 'Error opening specified xmol file for input'

  ! set OM_length to use
  select case (trim(style))
  case ('tip4p/2005')
     OM_length = 0.1546_dp
  case ('tip4p/Ice')
     OM_length = 0.1577_dp
  case ('tip4p')
     OM_length = 0.1500_dp
  case default
     write(*,*)'Unknown TIP4P variant. Select either:'
     write(*,*)' tip4p/2005'
     write(*,*)' tip4p/Ice'
     write(*,*)' tip4p'
     stop
  end select


  ! read the number of atoms - this must be the total number of O plus H
  ! sites - does not include the M sites.
  read(25,*,iostat=ierr) N
  if (ierr/=0) then
     stop 'Input xmolfile does not contain integer number of atoms on line 1'
  endif

  ! sometimes there's a character before the cell matrix
  read(25,*,iostat=ierr) dumchar, dumh           ! blank comment line in xmol format
  if (ierr/=0) then
     rewind(25)
     read(25,*)
     read(25,*,iostat=ierr) dumh
     if (ierr/=0) then
        stop 'Could not extract matrix of cell vectors from line 2 of xmolfile'
     end if
  end if
  
  ! specify cell dimensions here - NB orthorhombic only for now
  cell(1) = dumh(1,1)
  cell(2) = dumh(2,2)
  cell(3) = dumh(3,3)

  ! write header of config file
  write(*,'(a38)')'tip4p configuration - generated from config.xmol'
  write(*,'(2i10)')0,3
  write(*,'(3f20.10)')cell(1),0.0_dp,0.0_dp
  write(*,'(3f20.10)')0.0_dp,cell(2),0.0_dp
  write(*,'(3f20.10)')0.0_dp,0.0_dp,cell(3)

  ! initialise counters
  l = 1
  foundcount = 0     ! total number of molecules processed

  do i = 1,N

     rewind(25)      ! back to start of xmol file and 
     read(25,*) N    ! skip past header.
     read(25,*)

     found_this_pass = 0  ! found no oxygens yet on this pass

     do k = 1,N
        read(25,*,iostat=ierr)elemid,vec1
        if (ierr/=0) stop 'Unexpectedly reached end of xmolfile'
        if ((elemid/='H').and.(elemID/='O')) then
           write(*,*)'Error - found element of type '//elemid//' in xmol file.'
           write(*,*)'Input file must contain only oxygen and hydrogen.'
        end if
        if (elemid == 'O') then  ! found an oxygen on this pass
           found_this_pass = found_this_pass + 1
           !write(*,*)'Found ',found_this_pass,' oxygens in pass ',i
        end if
        if (found_this_pass>foundcount ) then
           foundcount = foundcount + 1  ! this is a new oxygen
           !write(*,*)'Found a new oxygen'
           exit
        end if
     end do

     rewind(25)
     read(25,*) N
     read(25,*)

     secondbond = .false.  ! no O-H bonds as yet

     do j = 1,N

        read(25,*)elemid,vec4
        if (elemid /= 'H') cycle  ! only interested in bonds to H atoms 

        oh1 = vec4 - vec1  ! compute separation vector using PBCs
        do idim = 1,3
           oh1(idim) = oh1(idim) - cell(idim)*anint(oh1(idim)/cell(idim))
        end do

        length = sqrt(dot_product(oh1,oh1))  ! length of this vector


        ! if less than 1.1 then add place H at specified length along
        ! that vector.
        if (length< 1.1_dp) then          
           if (secondbond) then
              vec3 = vec1 + oh1
              exit
           else
              vec2 = vec1 + oh1
              secondbond=.true.
           end if
        end if

     end do

     oh1 = vec2 - vec1
     oh2 = vec3 - vec1
     
     mag_oh1 = sqrt(dot_product(oh1,oh1))
     mag_oh2 = sqrt(dot_product(oh2,oh2))   

     ! and now adjust the two OH - lengths   
     mag_oh1 = sqrt(dot_product(oh1,oh1))
     mag_oh2 = sqrt(dot_product(oh2,oh2))   
     oh1 = oh1*OH_length/mag_oh1
     oh2 = oh2*OH_length/mag_oh2


     ! now need to find two perpendicular
     ! unit vectors in the plane of the molecule
     unit1 = oh1
     vec4  = cross_product(unit1,oh2)
     
     ! the second unit vector in the plane
     ! is simultaneously perpendicular to vec4 and unit1
     ! and hence is the cross-product of these two
     unit2 = cross_product(vec4,unit1)
     
     ! unitise
     unit2 = unit2/sqrt(dot_product(unit2,unit2))
     unit1 = unit1/sqrt(dot_product(unit1,unit1))     

     oh2 = unit2*sin(HOH_angle) + unit1*cos(HOH_angle)
     oh2 = oh2*OH_length

     vec2 = vec1 + oh1
     vec3 = vec1 + oh2

     ! check
     !costheta = dot_product(oh1,oh2)/(OH_length**2)
     !theta = acos(costheta)
     !write(*,*)'Angle H-O-H is ',theta*180.0_dp/Pi

     ! found both bonds - now add msite and print
     call m_position(vec1,vec2,vec3,vec4)

     ! print the molecule to the config file
     write(*,'(a7,5x,i6,5x,i6)')'O_tip4p',l,8
     write(*,'(3f20.10)')vec1
     write(*,'(a7,5x,i6,5x,i6)')'H_tip4p',l+1,1
     write(*,'(3f20.10)')vec2   
     write(*,'(a7,5x,i6,5x,i6)')'H_tip4p',l+2,1
     write(*,'(3f20.10)')vec3
     write(*,'(a7,5x,i6,5x,i6)')'M_tip4p',l+3,0
     write(*,'(3f20.10)')vec4     
     l = l + 4  ! used 4 atoms here

     if ( foundcount == N/3) exit  ! found all molecules

  end do

  close(25)

contains

  subroutine m_position(vec_o,vec_h1,vec_h2,vec_m)

    implicit none

    real(kind=dp),dimension(3),intent(in)  :: vec_o
    real(kind=dp),dimension(3),intent(in)  :: vec_h1
    real(kind=dp),dimension(3),intent(in)  :: vec_h2
    real(kind=dp),dimension(3),intent(out) :: vec_m
    real(kind=dp),dimension(3) :: res

    ! compute the vector along the bisector
    res = (vec_h1 - vec_o) + (vec_h2 - vec_o)

    ! unitise it
    res = res / sqrt(dot_product(res,res))

    ! position of the massless charge site is therefore..
    vec_m = vec_o + res*OM_Length


    return

  end subroutine m_position

  function cross_product(a,b)

    implicit none
    real(kind=dp),dimension(3),intent(in) :: a,b
    real(kind=dp),dimension(3) :: cross_product


    cross_product(1) = a(2)*b(3) - a(3)*b(2)
    cross_product(2) = -(a(1)*b(3)-a(3)*b(1))
    cross_product(3) = a(1)*b(2)-b(1)*a(2)

  end function cross_product


end program tip4p
