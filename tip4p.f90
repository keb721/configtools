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
!   5. Outputs as DL_POLY/GROMACS/LAMMPS input (LAMMPS with or without M-site)       !
!------------------------------------------------------------------------------------!
! version 2.0 - D.Quigley ( University of Warwick )                                  !   
!------------------------------------------------------------------------------------!
program tip4p

  implicit none

  integer,parameter :: dp = kind(1.0d0)      ! use double precision

  ! Pi and 1/3
  real(kind=dp),Parameter :: Pi = 3.141592653589793238462643383279502884197_dp,third = 1.0_dp/3

  ! integers/counters
  integer :: i,l,foundcount,found_this_pass,k,idim,N,j,mol

  ! vectors required for doing the geometry
  real(kind=dp),dimension(3)   :: vec1,vec2,vec3,vec4,cell,oh1,oh2
  real(kind=dp),dimension(3)   :: unit1,unit2

  ! temporary matrix of cell vectors
  real(kind=dp),dimension(3,3) :: dumh

  ! Distance OH common to all TIP4P variants 
  real(kind=dp),parameter :: OH_length = 0.9572_dp

  ! H - O - H bond angle required in radians
  real(kind=dp),parameter :: HOH_angle = 104.52_dp*Pi/180.0_dp
  
  ! Particle masses common to all TIP4P variants
  real(kind=dp),parameter :: O_mass = 15.9994,H_mass = 1.008

  ! Distance OM to use for the current TIP4P variant
  real(kind=dp) :: OM_length = 0.15_dp 
  
  ! Charges to use for the current TIP4P variant
  real(kind=dp) :: qO = -1.040_dp,qH = 0.520_dp,qM=0.0_dp 

  character(1) :: elemid   ! chemical symbol for the current atom
  character(1) :: dumchar  ! dummy character for cell comment line

  logical :: secondbond    ! is this the second O-H bond for the current H

  real(kind=dp) :: length          ! current bond length
  real(kind=dp) :: theta,costheta  ! current bond angle and cosine thereof
  real(kind=dp) :: mag_oh1,mag_oh2 ! length of O-H bonds before correction

  integer                        :: iarg,idata,ierr=0
  integer                        :: num_args
  character(30), dimension(0:10) :: command_line
  character(30)                  :: filename,style,output

  integer,parameter :: xml = 25 ! unit number for input xmol file  
  integer,parameter :: gro = 35 ! unit number for gromacs topology file (if used)
  
  
  ! check that there are two arguments.
  num_args = iargc()

  if (num_args/=3) then
     write(0,*)
     write(0,*) '                         T I P 4 P                          '
     write(0,*)
     write(0,*) '    Usage: tip4p <xmolfile> <tip4p_style> <output_style>    '
     write(0,*)
     write(0,*) '           D. Quigley - University of Warwick               '
     write(0,*)
     stop
  end if

  ! get command line arguments
  do iarg = 1, num_args
     call getarg(iarg,command_line(iarg))
  end do

  filename = trim(command_line(1))
  style    = trim(command_line(2))
  output   = trim(command_line(3))

  ! open specified input file, bail if it's not there
  open(unit=xml,file=trim(filename),status='old',iostat=ierr)
  if (ierr/=0) stop 'Error opening specified xmol file for input'

  ! set OM_length to use
  select case (trim(style))
  case ('tip4p/2005')
     OM_length =  0.1546_dp
     qO        = -1.1128_dp
     qH        =  0.5564_dp
  case ('tip4p/Ice')
     OM_length =  0.1577_dp
     qO        = -1.1794_dp
     qH        =  0.5897_dp
  case ('tip4p/long')
     OM_length =  0.1250_dp
     qO        = -1.0484_dp
     qH        =  0.5242_dp
  case ('tip4p')
     OM_length =  0.1500_dp
     qO        = -1.040_dp
     qH        =  0.520_dp
  case default
     write(*,*)'Unknown TIP4P variant. Select either:'
     write(*,*)' tip4p/2005'
     write(*,*)' tip4p/Ice'
     write(*,*)' tip4p/long [recommended parameters for use with long-range Coulombic solver]'
     write(*,*)' tip4p'
     stop
  end select

  ! read the number of atoms - this must be the total number of O plus H
  ! sites - does not include the M sites.
  read(xml,*,iostat=ierr)N
  if (ierr/=0) then
     stop 'Input xmolfile does not contain integer number of atoms on line 1'
  end if

  ! sometimes there's a character before the cell matrix
  read(xml,*,iostat=ierr)dumchar,dumh           ! blank comment line in xmol format
  if (ierr/=0) then
     rewind(xml)
     read(xml,*)
     read(xml,*,iostat=ierr)dumh
     if (ierr/=0) then
        stop 'Could not extract matrix of cell vectors from line 2 of xmolfile'
     end if
  end if
  
  ! specify cell dimensions here - NB orthorhombic only for now
  cell(1) = dumh(1,1)
  cell(2) = dumh(2,2)
  cell(3) = dumh(3,3)

  ! set output style to use
  select case (trim(output))
  case ('dlpoly')
     call write_dlpoly_header(style, cell)
  case ('lammps')
     call write_lammps_header(cell,style,N,O_mass,H_mass,.false.)
  case ('lammpsm')
     qM = qO
     qO = 0.0_dp
     call write_lammps_header(cell,style,int(4*N*third),O_mass,H_mass,.true.)
  case ('gromacs')
     call write_gromacs_header(style,int(4*N*third),gro)
  case default
     write(*,*)'Unknown output type. Select either:'
     write(*,*)' dlpoly'
     write(*,*)' lammpsm [LAMMPS output with explicit M-site]'
     write(*,*)' lammps  [LAMMPS output containing just O/H]'
     write(*,*)' gromacs [gro file written to command line, topology file output to tip4p_gromacs.top]'
     stop
  end select

  ! initialise counters
  l          = 1
  mol        = 1
  foundcount = 0     ! total number of molecules processed

  do i = 1, N

     rewind(xml)      ! back to start of xmol file and 
     read(xml,*)N    ! skip past header.
     read(xml,*)

     found_this_pass = 0  ! found no oxygens yet on this pass

     do k = 1,N
        read(xml,*,iostat=ierr)elemid,vec1
        if (ierr/=0) stop 'Unexpectedly reached end of xmolfile'
        if ((elemid/='H').and.(elemID/='O')) then
           write(*,*)'Error - found element of type '//elemid//' in xmol file.'
           write(*,*)'Input file must contain only oxygen and hydrogen.'
        end if
        if (elemid=='O') then  ! found an oxygen on this pass
           found_this_pass = found_this_pass + 1
           !write(*,*)'Found ',found_this_pass,' oxygens in pass ',i
        end if
        if (found_this_pass>foundcount) then
           foundcount = foundcount + 1  ! this is a new oxygen
           !write(*,*)'Found a new oxygen'
           exit
        end if
     end do

     rewind(xml)
     read(xml,*)N
     read(xml,*)

     secondbond = .false.  ! no O-H bonds as yet

     do j = 1,N

        read(xml,*)elemid,vec4
        if (elemid/='H') cycle  ! only interested in bonds to H atoms 

        oh1 = vec4 - vec1  ! compute separation vector using PBCs
        do idim = 1, 3
           oh1(idim) = oh1(idim) - cell(idim)*anint(oh1(idim)/cell(idim))
        end do

        length = sqrt(dot_product(oh1,oh1))  ! length of this vector


        ! if less than 1.1 then add place H at specified length along
        ! that vector.
        if (length<1.1_dp) then          
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

     if (trim(output)=='dlpoly') then
        call m_position(vec1,vec2,vec3,vec4)
        call write_dlpoly_config(l,vec1,vec2,vec3,vec4)
     else if (trim(output)=='lammps') then
        call write_lammps_config(l,mol,vec1,vec2,vec3,vec4,qO,qH,qM,.false.)
     else if (trim(output)=='lammpsm') then
        call m_position(vec1,vec2,vec3,vec4)
        call write_lammps_config(l,mol,vec1,vec2,vec3,vec4,qO,qH,qM,.true.)
     else if (trim(output)=='gromacs') then
        call m_position(vec1,vec2,vec3,vec4)
        call write_gromacs_config(l,mol,vec1,vec2,vec3,vec4,qO,qH,O_mass,H_mass,gro)
     end if
       
     if (foundcount==N/3) exit  ! found all molecules

  end do

  close(xml)
  
  if (trim(output)=='lammps') then
     call write_lammps_footer(N)
  else if (trim(output)=='lammpsm') then
     call write_lammps_footer_m(int(4*N*third))
  else if (trim(output)=='gromacs') then
     call write_gromacs_footer(int(4*N*third),cell,gro)
  end if

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
    res = res/sqrt(dot_product(res,res))

    ! position of the massless charge site is therefore..
    vec_m = vec_o + res*OM_Length


    return

  end subroutine m_position

  function cross_product(a,b)

    implicit none
    real(kind=dp),dimension(3),intent(in) :: a,b
    real(kind=dp),dimension(3) :: cross_product

    cross_product(1) =   a(2)*b(3) - a(3)*b(2)
    cross_product(2) = -(a(1)*b(3) - a(3)*b(1))
    cross_product(3) =   a(1)*b(2) - b(1)*a(2)

  end function cross_product
  
  subroutine write_lammps_header(cell,style,N,O_mass,H_mass,Ms)

    implicit none
    real(kind=dp),dimension(3),intent(in) :: cell
    real(kind=dp),intent(in)              :: O_mass,H_mass
    character(30),intent(in)              :: style
    integer,intent(in)                    :: N
    logical,intent(in)                    :: Ms
    integer                               :: mols
    
     ! write header of config file
     write(*,'(a38)')'# '//trim(style)//' configuration'

     write(*,*)
     write(*,'(f12.8,a1,f12.8,a8)')0.0,' ',cell(1),' xlo xhi' 
     write(*,'(f12.8,a1,f12.8,a8)')0.0,' ',cell(2),' ylo yhi'  
     write(*,'(f12.8,a1,f12.8,a8)')0.0,' ',cell(3),' zlo zhi'  
     write(*,*)

     write(*,*)
     write(*,'(i4,a6)')N,' atoms'

     if (Ms) then
        mols = 0.25*N
        write(*,'(i4,a6)')3*mols,' bonds'  
        write(*,'(i4,a7)')3*mols,' angles'
        write(*,*)
        write(*,'(i1,a11)')3,' atom types'
        write(*,'(i1,a11)')2,' bond types' 
        write(*,'(i1,a12)')2,' angle types'

     else    
        mols = third*N
        write(*,'(i4,a6)')2*mols,' bonds'
        write(*,'(i4,a7)')mols,' angles'
        write(*,*)
        write(*,'(i1,a11)')2,' atom types'
        write(*,'(i1,a11)')1,' bond types' 
        write(*,'(i1,a12)')1,' angle types'
     end if

     write(*,*) 
     ! define atoms with their masses
     write(*,*) 
     write(*,*)'Masses'
     write(*,*)
     write(*,'(i1,a1,f12.8,a4)')1,' ', O_mass, ' # O'
     write(*,'(i1,a1,f12.8,a4)')2,' ', H_mass, ' # H'
     if (Ms) write(*,'(i1,a1,f12.8,a4)')3,' ',0.00000001, ' # M'   ! LAMMPS will not accept a truly massless M-site, therefore make mass very small

     write(*,*)

     ! write atomic positions
     write(*,*)
     write(*,*)'Atoms'
     write(*,*)'# id mol type  charge         x            y            z ' 
     
    return

  end subroutine write_lammps_header

  subroutine write_lammps_config(l,mol,vec1,vec2,vec3,vec4,qO,qH,qM,Ms)

     implicit none
     integer,intent(inout)                 :: l,mol
     real(kind=dp),dimension(3),intent(in) :: vec1,vec2,vec3,vec4     
     real(kind=dp),intent(in)              :: qO,qH,qM
     logical,intent(in)                    :: Ms
     character(51)                         :: form = '(i5,a1,i5,a2,i1,a2,f9.6,a1,f12.8,a1,f12.8,a1,f12.8)'
     
     write(*,form)l,' ',mol,'  ',1,'  ',qO,' ',vec1(1),' ',vec1(2),' ',vec1(3)
     write(*,form)l+1,' ',mol,'  ',2,'  ',qH,' ',vec2(1),' ',vec2(2),' ',vec2(3)
     write(*,form)l+2,' ',mol,'  ',2,'  ',qH,' ',vec3(1),' ',vec3(2),' ',vec3(3)
     l = l + 3
     if (Ms) then
        write(*,form)l,' ',mol,'  ',3,'  ',qM,' ',vec4(1),' ',vec4(2),' ',vec4(3)
        l = l + 1 
     end if

     mol = mol + 1

    return

  end subroutine write_lammps_config

  subroutine write_lammps_footer(N)

    implicit none
    integer,intent(in) :: N
    integer            :: i, atn

    write(*,*)
    write(*,*)
    
    ! Write bonds information
    write(*,*)'Bonds'
    write(*,*)'# bondid type atm1 atm2'

    i = 1
    do atn = 1, N, 3
       write(*,'(i8,a2,i1,a2,i4,a1,i4)')i,'  ',1,'  ',atn,' ', atn+1
       write(*,'(i8,a2,i1,a2,i4,a1,i4)')i+1,'  ',1,'  ',atn,' ', atn+2
       i = i + 2
    end do
       
    write(*,*)
    write(*,*)
    
    ! Write angles information
    write(*,*)'Angles'
    write(*,*)'# anglid type atm1 atm2 atm3'
    i = 1
    do atn = 1, N, 3
       write(*,'(i8,a2,i1,a2,i4,a1,i4,a1,i4)')i,'  ',1,'  ',atn+1,' ',atn,' ',atn+2 
       i = i + 1
    end do
    
    return

  end subroutine write_lammps_footer

  subroutine write_lammps_footer_m(N)

    implicit none
    integer,intent(in) :: N
    integer            :: i, atn

    write(*,*)
    write(*,*)
    
    ! Write bonds information
    write(*,*)'Bonds'
    write(*,*)'# bondid type atm1 atm2'

    i = 1
    do atn = 1, N, 4
       write(*,'(i8,a2,i1,a2,i4,a1,i4)')i,'  ',1,'  ',atn,' ', atn+1
       write(*,'(i8,a2,i1,a2,i4,a1,i4)')i+1,'  ',1,'  ',atn,' ', atn+2
       write(*,'(i8,a2,i1,a2,i4,a1,i4)')i+2,'  ',2,'  ',atn,' ', atn+3
       i = i + 3
    end do
       
    write(*,*)
    write(*,*)
    
    ! Write angles information
    write(*,*)'Angles'
    write(*,*)'# anglid type atm1 atm2 atm3'
    i = 1
    do atn = 1, N, 4
       write(*,'(i8,a2,i1,a2,i4,a1,i4,a1,i4)')i,'  ',1,'  ',atn+1,' ',atn,' ',atn+2 
       write(*,'(i8,a2,i1,a2,i4,a1,i4,a1,i4)')i+1,'  ',2,'  ',atn+1,' ',atn,' ',atn+3
       write(*,'(i8,a2,i1,a2,i4,a1,i4,a1,i4)')i+2,'  ',2,'  ',atn+2,' ',atn,' ',atn+3 
       i = i + 3
    end do

    return 

  end subroutine write_lammps_footer_m

  subroutine write_gromacs_header(style,N,gro)

    implicit none
    character(30),intent(in)              :: style
    integer,intent(in)                    :: N,gro
    integer                               :: ierr

    write(*,'(a38)')trim(style)//' configuration'
    write(*,'(i4)')N

    ! ================ !
    ! WRITING TOP FILE ! 
    ! ================ !
    
    open(gro,file='tip4p_gromacs.top',status='new',iostat=ierr)
    if (ierr/=0) stop 'Error opening `tip4p_groamcs.top` file for output'

    write(gro,*)';'
    write(gro,*)'; Topology file written by refine_bonds.f90 from atomic hash'
    write(gro,*)';'
    write(gro,*)'; Force-field files to be included'
    write(gro,*)'#include '
    write(gro,*)
    write(gro,*)
    write(gro,*)'[ moleculetype ]'
    write(gro,*)'; name nrexcl'
    write(gro,*)'SOL 1'
    write(gro,*)
    write(gro,*)
    write(gro,*)'[ atoms ]'
    write(gro,*)'; counter atomtype residuenum residuename atomname chargegroup charge mass'

    return

  end subroutine write_gromacs_header


  subroutine write_gromacs_config(l,mol,vec1,vec2,vec3,vec4,qO,qH,O_mass,H_mass,gro)

    implicit none
    integer,intent(inout)                 :: l,mol
    real(kind=dp),dimension(3),intent(in) :: vec1,vec2,vec3,vec4     
    real(kind=dp),intent(in)              :: qO,qH,O_mass,H_mass
    integer,intent(in)                    :: gro
    real(kind=dp)                         :: cf = 0.1_dp ! Gromacs output in nm, not A as DL_POLY and LAMMPS
    character(52)                         :: formt = '(i4,a1,a2,a1,i4,a1,a5,a1,a5,a1,i4,a1,f12.8,a1,f12.8)'

    write(*,'(i5,2a5,i5,3F8.3)')mol,'SOL','OW',l,vec1(1)*cf,vec1(2)*cf,vec1(3)*cf
    write(*,'(i5,2a5,i5,3F8.3)')mol,'SOL','HW1',l+1,vec2(1)*cf,vec2(2)*cf,vec2(3)*cf
    write(*,'(i5,2a5,i5,3F8.3)')mol,'SOL','HW2',l+2,vec3(1)*cf,vec3(2)*cf,vec3(3)*cf
    write(*,'(i5,2a5,i5,3F8.3)')mol,'SOL','MW3',l+3,vec4(1)*cf,vec4(2)*cf,vec4(3)*cf

    write(gro,formt)l,' ','O',' ',mol,' ','SOL',' ','OW',' ',l,' ',0.0,' ',O_mass
    write(gro,formt)l+1,' ','H',' ',mol,' ','SOL',' ','HW1',' ',l+1,' ',qH,' ',H_mass
    write(gro,formt)l+2,' ','H',' ',mol,' ','SOL',' ','HW2',' ',l+2,' ',qH,' ',H_mass
    write(gro,formt)l+3,' ','M',' ',mol,' ','SOL',' ','MW',' ',l+3,' ',qO,' ',0.0    

    l   = l + 4
    mol = mol + 1

    return

  end subroutine write_gromacs_config


  subroutine write_gromacs_footer(N,cell,gro)
    
    implicit none
    real(kind=dp),dimension(3),intent(in) :: cell
    integer,intent(in)                    :: N,gro
    real(kind=dp)                         :: cf = 0.1_dp ! Gromacs output in nm, not A as DL_POLY and LAMMPS
    integer                               :: i, j
    
    write(*,'(f12.8,a1,f12.8,a1,f12.8)') cell(1)*cf,' ',cell(2)*cf,' ',cell(3)*cf

    ! End of .gro file

    write(gro,*)

    write(gro,*)'[ bonds ]'
    write(gro,*)'; atom1 atom2'
    
    do i = 1, N, 4
       do j = 1, 3
          write(gro,'(i4,a1,i4)')i,' ',i+j
       end do
    end do
    
    write(gro,*)
    
    ! Write angles information
    write(gro,*)'[ angles ]'
    write(gro,*)'; atom1 atom2 atom3'
    
    do i = 1, N, 4
       write(gro,'(i4,a1,i4,a1,i4)')i+1,' ',i,' ',i+2
       write(gro,'(i4,a1,i4,a1,i4)')i+1,' ',i,' ',i+3
       write(gro,'(i4,a1,i4,a1,i4)')i+2,' ',i,' ',i+3
    end do

    close(gro)

    return
    
  end subroutine write_gromacs_footer


  subroutine write_dlpoly_header(style, cell)

    implicit none
    real(kind=dp),dimension(3),intent(in)   :: cell
    character(30),intent(in)                :: style
    
    
     ! write header of config file
     write(*,'(a38)')trim(style)//' configuration'
     write(*,'(2i10)')0,3
     write(*,'(3f20.10)')cell(1),0.0_dp,0.0_dp
     write(*,'(3f20.10)')0.0_dp,cell(2),0.0_dp
     write(*,'(3f20.10)')0.0_dp,0.0_dp,cell(3)

    return

  end subroutine write_dlpoly_header

  subroutine write_dlpoly_config(l,vec1,vec2,vec3,vec4)

    implicit none
    integer,intent(inout)                 :: l
    real(kind=dp),dimension(3),intent(in) :: vec1,vec2,vec3,vec4
    
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

    return 

  end subroutine write_dlpoly_config




end program tip4p
