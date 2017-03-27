program zif_cif2gin
 use iso_fortran_env
 implicit none
! locals
 integer             :: i,j,k,ierr,nn
! variables
 integer             :: num_args
 integer             :: n_atoms = 0
 real                :: cell_0(1:6) = 0.0
 real                :: rv(3,3),vr(3,3)
 character(len=100)  :: line
 character(len=20)   :: spam
 character(len=80)   :: string_stop_head= "_atom_site_charge" !"_atom_site_charge"
 character(len=100)  :: CIFFilename=" "
 character(len=100)  :: filename=" "
 ! allocatable
 real,allocatable    :: xcrystal(:,:),xcartes(:,:),charge(:)
 real,allocatable    :: DistanceMatrix(:,:)
 character(len=100), dimension(:), allocatable :: args
 logical,dimension(:,:),allocatable            :: adj
 character(len=4),allocatable                  :: label(:)
 num_args = command_argument_count()
 allocate(args(num_args))
 do i = 1, num_args
  call get_command_argument(i,args(i))
 end do
 write(6,'(a,1x,i2)')'Arguments:',command_argument_count()
 write(6,'(a)')( args(i),i=1,num_args)
 if(num_args==0) then
  call print_help()
  stop
 end if
 do i=1,num_args
  select case(args(i))
   case ('-h','--help')
    call print_help()
    stop
   case ('-c','--cif')
    CIFFilename=args(i+1)
    filename=CIFFilename(1:Clen_trim(CIFFilename)-4)
    write(6,'(a)') filename
  end select
 end do
 beta = 1.0/(k_B*temperature)
 open(100,file=CIFfilename,status='old',iostat=ierr)
 if(ierr/=0) stop 'Error opening CIF file'
 read_cif: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit read_cif
  if(line(1:14)=="_cell_length_a")then
   read(line,*)spam,cell_0(1)
   cycle read_cif
  end if
  if(line(1:14)=="_cell_length_b")then
   read(line,*)spam,cell_0(2)
   cycle read_cif
  end if
  if(line(1:14)=="_cell_length_c")then
   read(line,*)spam,cell_0(3)
   cycle read_cif
  end if
  if(line(1:17)=="_cell_angle_alpha")then
   read(line,*)spam,cell_0(4)
   cycle read_cif
  end if
  if(line(1:16)=="_cell_angle_beta")then
   read(line,*)spam,cell_0(5)
   cycle read_cif
  end if
  if(line(1:17)=="_cell_angle_gamma")then
   read(line,*)spam,cell_0(6)
   cycle read_cif
  end if
  if(line(1:)==string_stop_head) exit read_cif
 end do read_cif
 call cell(rv,vr,cell_0)
 n_atoms=0
 read_natoms: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit read_natoms
  n_atoms=n_atoms+1
 end do read_natoms
!
 allocate(xcrystal(3,n_atoms),xcartes(3,n_atoms))
 allocate(charge(n_atoms))
 allocate(DistanceMatrix(n_atoms,n_atoms))
 allocate(label(n_atoms))
 allocate(adj(n_atoms,n_atoms))
!
 rewind(100)
 write(6,*)'Atoms:',n_atoms
 do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit
  if(line(1:)==string_stop_head) exit
 end do
 i=0
 read_atoms: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit
  i=i+1
  read(line,*)label(i),( xcrystal(j,i),j=1,3),charge(i)
  do j=1,3
   xcrystal(j,i)=MOD(xcrystal(j,i)+100.0,1.0)
   xcartes(j,i) =rv(j,1)*xcrystal(1,i) + rv(j,2)*xcrystal(2,i)+rv(j,3)*xcrystal(3,i)
  end do
  write(6,'(a4,1x,3(f14.12,1x),3(f14.7,1x))',ADVANCE='yes') &
   label(i),( xcrystal(j,i),j=1,3),( xcartes(j,i),j=1,3)
 end do read_atoms
 call make_dist_matrix(n_atoms,cell_0,rv,vr,xcrystal,DistanceMatrix)
 write(6,*)'=========='
 close(100)
! Analysis:
 call write_graph()
!
 deallocate(xcrystal)
 deallocate(xcartes)
 deallocate(DistanceMatrix)
 deallocate(charge)
 deallocate(label)
 deallocate(newlabel)
 contains
 subroutine write_graph()
  implicit none
  integer  i
 end subroutine write_graph
 PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
 CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
 INTEGER :: i
 Clen = LEN(s)
 i = LEN_TRIM(s)
 IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen
 PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
 CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
 INTEGER :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
 i = LEN_TRIM(s) ; Clen_trim = i
 IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim
!
 SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real, parameter   :: pi = ACOS(-1.0)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: DEGTORAD
 DEGTORAD=pi/180.0
 IF(cell_0(4) == 90.0) THEN
   cosa = 0.0
 ELSE
   ALP=cell_0(4)*degtorad
   COSA=cos(ALP)
 ENDIF
 IF(cell_0(5) == 90.0) THEN
   cosb = 0.0
 ELSE
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 ENDIF
 IF(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 ELSE
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 ENDIF
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
 call inverse(rv,vr,3)
 WRITE(6,'(a)') 'Cell:'
 WRITE(6,'(6F14.7)')( cell_0(j), j=1,6 )
 WRITE(6,'(a)')'Linear Transformation Operator:'
 DO i=1,3
  WRITE(6,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 WRITE(6,'(a)')'Inverse Linear Transformation Operator:'
 DO i=1,3
  WRITE(6,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 RETURN
 END SUBROUTINE cell
!
 SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  REAL               :: radtodeg
  REAL, PARAMETER    :: pi=ACOS(-1.0) 
  radtodeg=180.0/PI
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    enddo
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
  RETURN
 END SUBROUTINE uncell
!
 SUBROUTINE inverse(a,c,n)
 implicit none
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 L=0.0
 U=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  L(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 RETURN
 END SUBROUTINE inverse
!
subroutine make_dist_matrix(n,cell_0,rv,vr,x,dist_matrix)
 implicit none
 integer,intent(in) :: n
 real,intent(in)    :: cell_0(6),rv(3,3),vr(3,3),x(3,n)
 real,intent(out)   :: dist_matrix(n,n)
 integer            :: i,j,k
 real               :: r1(3),r2(3),s
 DO i=1,n
    dist_matrix(i,i)=0.0
    DO j=i+1,n
       forall ( k=1:3 )
        r1(k)=x(k,i)
        r2(k)=x(k,j)
       end forall
       call make_distances(cell_0,r1,r2,rv,s)
       dist_matrix(i,j)=s
       dist_matrix(j,i)=dist_matrix(i,j)
    END DO
 END DO
 return
end subroutine make_dist_matrix
!
 SUBROUTINE make_distances(cell_0,r2,r1,rv,dist)
 IMPLICIT NONE
 REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist
 REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
 INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
 REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
   enddo
  enddo
  dist=MINVAL(d_image)
  RETURN
 END SUBROUTINE
!
 REAL FUNCTION DISTANCE(atom,ouratom,rv)
  IMPLICIT NONE
  INTEGER :: j
  REAL :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
  REAL :: rv(3,3)
  FORALL ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
  END FORALL
  DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
 END FUNCTION
end program zif_cif2gin

