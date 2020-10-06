module functions
  implicit none
contains

function inverse(A)
  implicit none
  integer(1) :: i,j,N
  real(8) :: inverse(3,3)
  real(8) :: A(:,:)
  real(8) :: B(3,3),IA(3,3)

  N=3
  forall(i=1:N,j=1:N,i==j) IA(i,j)=1.D0
  forall(i=1:N,j=1:N,i/=j) IA(i,j)=0.D0
  B=A
  call upper(B,IA,N)
  call lower(B,IA,N)
  forall(i=1:N) IA(i,:)=IA(i,:)/B(i,i)
  inverse=IA
  return
end function inverse

subroutine upper(M,S,N)
  implicit none
  integer(1) :: N
  real(8) :: M(N,N)
  real(8) :: S(N,N)
  integer(1) :: I,J
  real(8) :: E

  do I=1,N-1
    do J=I+1,N
      E=M(J,I)/M(I,I)
      M(J,I:N)=M(J,I:N)-M(I,I:N)*E
      S(J,:)=S(J,:)-S(I,:)*E
    end do
  end do
  return
end subroutine upper

subroutine lower(M,S,N)
  implicit none
  integer(1) :: N
  real(8) :: M(N,N)
  real(8) :: S(N,N)
  integer(1) :: I,J
  real(8) :: E

  do I=N,2,-1
    do J=I-1,1,-1
      E=M(J,I)/M(I,I)
      M(J,I:N)=M(J,I:N)-M(I,I:N)*E
      S(J,:)=S(J,:)-S(I,:)*E
    end do
  end do
  return
end subroutine lower

function det(A)
  implicit none
  integer(1) :: det,A(3,3)

  det=    A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
  det=det-A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  det=det+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  return
end function det

function array(start,step,finish)
  implicit none
  real(8),allocatable :: array(:)
  real(8) :: start,step,finish
  integer(2) :: i,n

  n=nint((finish-start)/step+1)
  allocate(array(n))
  forall(i=1:n)
    array(i)=start+step*(i-1)
  end forall
  return
end function array

function cross(A,B)
  implicit none
  real(8) :: cross(3)
  real(8) :: A(3),B(3)

  cross(1)=A(2)*B(3)-A(3)*B(2)
  cross(2)=A(3)*B(1)-A(1)*B(3)
  cross(3)=A(1)*B(2)-A(2)*B(1)
  return
end function cross

end module functions

module symmetry
  implicit none
contains

subroutine rotation(rot,a,x,natom,symprec)
  use functions
  implicit none
  real(8) :: a(3,3),x(:,:),symprec
  integer(2) :: natom(:),satom
  real(8) :: G(3,3),dG(3,3),x_t(3),T(3),dx(3),c(3)
  integer(1) :: W11,W12,W13,W21,W22,W23,W31,W32,W33
  integer(1) :: W(3,3),rot_t(3,3,48)
  integer(2) :: i,j,k,m,n,m1,m2,n1,n2,num
  integer(1),allocatable :: rot(:,:,:)
  logical(1) :: discard

  G=matmul(transpose(a),a)
  n=0
  do W11=-1,1
    do W12=-1,1
      do W13=-1,1
        do W21=-1,1
          do W22=-1,1
            do W23=-1,1
              do W31=-1,1
                do W32=-1,1
                  do W33=-1,1
                    W(1,:)=(/ W11, W12, W13 /)
                    W(2,:)=(/ W21, W22, W23 /)
                    W(3,:)=(/ W31, W32, W33 /)
                    if ( abs(det(W)) == 1 ) then
                      dG=dabs(G-matmul(matmul(transpose(W),G),W))
                      if ( all(dG < symprec) ) then
                        n=n+1
                        rot_t(:,:,n)=W
                      end if
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  m=minloc(natom,1)
  n1=natom(m)
  m2=sum(natom(1:m))
  m1=m2-n1+1
  satom=sum(natom)
  num=0
  do i=1,n
    W=rot_t(:,:,i)
    do j=m1,m2
      T=x(:,j)-matmul(W,x(:,m1))
      k=1
      n1=1
      n2=natom(1)
      do m=1,satom
        x_t=matmul(W,x(:,m))+T
        discard=.true.
        do n=n1,n2
          dx=x_t-x(:,n)
          c=dabs(matmul(a,(dx-nint(dx))))
          if ( all(c < symprec) ) then
            discard=.false.
            exit
          end if
        end do
        if ( discard == .true. ) then
          exit
        end if        
        if ( m == sum(natom(1:k)) .and. m /= satom ) then
          k=k+1
          n1=n2+1
          n2=n2+natom(k)
        end if
      end do
      if ( discard == .false. ) then
        num=num+1
        rot_t(:,:,num)=W
        exit
      end if
    end do
  end do

  allocate(rot(3,3,num))
  rot=rot_t(:,:,1:num)
  return
end subroutine rotation

end module symmetry

module structure
  implicit none
contains

subroutine readpos(a,x,natom,poscar)
  use functions
  implicit none
  real(8) :: a(3,3)
  real(8),allocatable :: x(:,:)
  integer(1) :: ntype,io_error
  integer(2) :: satom,i
  integer(2),allocatable :: natom(:)
  real(8) :: scal
  character(len=1) :: coordinate
  character(len=*) :: poscar

  open(10,file=poscar)
  read(10,*)
  read(10,*) scal
  read(10,*) a
  a=scal*a
  read(10,*)
  ntype=1
  do while (io_error <= 0)
    allocate(natom(ntype))
    read(10,*,iostat=io_error) natom
    ntype=ntype+1
    deallocate(natom)
    backspace(10)
  end do
  ntype=ntype-2
  allocate(natom(ntype))
  backspace(10)
  read(10,*) natom
  satom=sum(natom)
  allocate(x(3,satom))
  read(10,*) coordinate
  read(10,*) x
  close(10)
  if ( coordinate == 'C' ) then
    x=matmul(inverse(a),x)
  end if
  return
end subroutine readpos

end module structure

module kpoints
  implicit none
contains

subroutine ibzkpt(ikpt,deg,nk,a,rot,gam,soc,symprec)
  use functions
  implicit none
  real(kind=8),parameter :: pi=3.141592653589793D0
  real(8),allocatable :: ikpt(:,:),akpt(:,:),ikpt_t(:,:)
  integer(1),allocatable :: deg(:),deg_t(:),sym(:,:,:)
  integer(1) :: nk(:),km,rot(:,:,:),inv(3,3),no
  logical(1) :: gam,soc
  real(8) :: a(3,3),b(3,3),dk(3),omega,okpt(3),dkpt(3)
  real(8) :: kx(nk(1)),ky(nk(2)),kz(nk(3)),symprec
  integer(2) :: ind,nkpts,i,j,k,n
  logical(1),allocatable :: occ(:)

  dk=nk
  kx=array(0.D0,1/dk(1),(dk(1)-1)/dk(1))
  ky=array(0.D0,1/dk(2),(dk(2)-1)/dk(2))
  kz=array(0.D0,1/dk(3),(dk(3)-1)/dk(3))

  km=(nk(1)+3)/2
  if ( nk(1) ) then
    if ( nk(1) > 1 ) kx(km:)=kx(km:)-1.D0
  else
    if ( gam == .true. ) then
      if ( nk(1) > 2 ) kx(km+1:)=kx(km+1:)-1.D0
    else
      kx=kx+1/(2*dk(1))
      kx(km:)=kx(km:)-1.D0
    end if
  end if

  km=(nk(2)+3)/2
  if ( nk(2) ) then
    if ( nk(2) > 1 ) ky(km:)=ky(km:)-1.D0
  else
    if ( gam == .true. ) then
      if ( nk(2) > 2 ) ky(km+1:)=ky(km+1:)-1.D0
    else
      ky=ky+1/(2*dk(2))
      ky(km:)=ky(km:)-1.D0
    end if
  end if

  km=(nk(3)+3)/2
  if ( nk(3) ) then
    if ( nk(3) > 1 ) kz(km:)=kz(km:)-1.D0
  else
    if ( gam == .true. ) then
      if ( nk(3) > 2 ) kz(km+1:)=kz(km+1:)-1.D0
    else
      kz=kz+1/(2*dk(3))
      kz(km:)=kz(km:)-1.D0
    end if
  end if

  nkpts=nk(1)*nk(2)*nk(3)
  allocate(akpt(3,nkpts))
  do i=0,nk(3)-1
    do j=0,nk(2)-1
        ind=i*(nk(1)*nk(2))+j*nk(1)
        akpt(1,ind+1:ind+nk(1))=kx
        akpt(2,ind+1:ind+nk(1))=spread(ky(j+1),1,nk(1))
    end do
    ind=i*(nk(1)*nk(2))
    akpt(3,ind+1:ind+nk(1)*nk(2))=spread(kz(i+1),1,nk(1)*nk(2))
  end do

  if ( soc == .true. ) then
    allocate(ikpt(3,nkpts))
    allocate(deg(nkpts))
    ikpt=akpt
    deg=1
    return
  end if

  omega=dot_product(cross(a(:,1),a(:,2)),a(:,3))
  b(:,1)=2*pi*cross(a(:,2),a(:,3))/omega
  b(:,2)=2*pi*cross(a(:,3),a(:,1))/omega
  b(:,3)=2*pi*cross(a(:,1),a(:,2))/omega

  inv(1,:)=(/ -1,  0,  0 /)
  inv(2,:)=(/  0, -1,  0 /)
  inv(3,:)=(/  0,  0, -1 /)
  no=size(rot,3)
  allocate(sym(3,3,2*no))
  do i=1,no
    sym(:,:,i)=transpose(rot(:,:,i))
    sym(:,:,no+i)=matmul(sym(:,:,i),inv)
  end do
  no=no*2

  allocate(ikpt_t(3,nkpts))
  allocate(deg_t(nkpts))
  allocate(occ(nkpts))
  occ=.true.
  deg_t=1

  n=0
  do i=1,nkpts
    if ( occ(i) == .false. ) cycle
    n=n+1
    ikpt_t(:,n)=akpt(:,i)
    do j=1,no
      okpt=matmul(sym(:,:,j),akpt(:,i))
      do k=i+1,nkpts
        if ( occ(k) == .true. ) then
          dkpt=akpt(:,k)-okpt
          dkpt=dabs(matmul(b,(dkpt-nint(dkpt))))
          if ( all(dkpt < symprec) ) then
            occ(k)=.false.
            deg_t(n)=deg_t(n)+1
          end if
        end if
      end do
    end do
  end do

  allocate(ikpt(3,n))
  allocate(deg(n))
  ikpt=ikpt_t(:,1:n)
  deg=deg_t(1:n)
  return
end subroutine ibzkpt

subroutine bandkpt(bkpt,hkpt,nk,ng)
  implicit none
  real(8),allocatable :: bkpt(:,:)
  real(8) :: hkpt(:,:,:),dkpt(3)
  integer(2) :: nk,nkpts,j
  integer(1) :: ng,i

  nkpts=nk*ng
  allocate(bkpt(3,nkpts))
  do i=1,ng
    dkpt=(hkpt(:,2,i)-hkpt(:,1,i))/(nk-1)
    do j=1,nk
      bkpt(:,(i-1)*nk+j)=hkpt(:,1,i)+dkpt*(j-1)
    end do
  end do
  return
end subroutine bandkpt

end module kpoints

program hsekpt
  use symmetry
  use structure
  use kpoints
  implicit none
  logical(1) :: alive,gam,soc
  integer(1) :: k(3),ng
  integer(2) :: nk,nkpts,i
  real(8) :: a(3,3)
  real(8),allocatable :: ikpt(:,:),bkpt(:,:)
  real(8),allocatable :: hkpt(:,:,:),x(:,:)
  integer(2),allocatable :: natom(:)
  integer(1),allocatable :: rot(:,:,:),deg(:)

  inquire(file='inkpt',exist=alive)
  if ( alive == .false. ) then
    open(11,file='inkpt')
    write(11,*) ' .false.  .false.  ! gamma  lsorbit'
    write(11,*)
    write(11,*) ' 5  5  1'
    write(11,*)
    write(11,*) ' 20  4'
    write(11,*)
    write(11,*) ' 0.000  0.000  0.000'
    write(11,*) '-0.333  0.667  0.000'
    write(11,*)
    write(11,*) '-0.333  0.667  0.000'
    write(11,*) ' 0.000  0.500  0.000'
    write(11,*)
    write(11,*) ' 0.000  0.500  0.000'
    write(11,*) ' 0.333  0.333  0.000'
    write(11,*)
    write(11,*) ' 0.333  0.333  0.000'
    write(11,*) ' 0.000  0.000  0.000'
    close(11)
    write(*,*) 'ERROR: inkpt file does not exist !'
    write(*,*) 'A sample file has been generated automatically,'
    write(*,*) 'Please modify it as your needed and try again.'
    stop
  end if

  open(12,file='inkpt')
  read(12,*) gam,soc
  read(12,*)
  read(12,*) k
  read(12,*)
  read(12,*) nk,ng
  allocate(hkpt(3,2,ng))
  do i=1,ng
    read(12,*)
    read(12,*) hkpt(:,1,i)
    read(12,*) hkpt(:,2,i)
  end do
  close(12)

  call readpos(a,x,natom,'POSCAR')
  call rotation(rot,a,x,natom,1D-4)
  call ibzkpt(ikpt,deg,k,a,rot,gam,soc,1D-3)
  call bandkpt(bkpt,hkpt,nk,ng) 

  nkpts=size(ikpt,2)
  open(13,file='IBZKPT')
  write(13,'(A)') 'Automatically generated mesh'
  write(13,'(I8)') nkpts
  write(13,'(A)') 'Reciprocal lattice'
  do i=1,nkpts
    write(13,'(3F20.14,I14)') ikpt(:,i),deg(i)
  end do
  close(13)

  open(14,file='KPOINTS')
  write(14,'(A)') 'IBZKPT and BANDKPT for HSE Calculations'
  write(14,'(I8)') nkpts+nk*ng
  write(14,'(A)') 'Reciprocal lattice'
  do i=1,nkpts
    write(14,'(3F20.14,I14)') ikpt(:,i),deg(i)
  end do
  do i=1,nk*ng
    write(14,'(3F20.14,I14)') bkpt(:,i),0
  end do
  close(14)
end program hsekpt
