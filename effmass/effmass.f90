module eig2band
  implicit none
contains

subroutine get_band(eigen,kpath,nkpts,nk)
  implicit none
  integer :: nkpts,nk,nbands,temp
  real(kind=8),allocatable :: eigen(:,:)
  real(kind=8),allocatable :: kpath(:)
  integer :: i,j

  open(12,file='EIGENVAL')
  read(12,*);read(12,*);read(12,*);read(12,*);read(12,*)
  read(12,*) temp,nkpts,nbands
  allocate(eigen(nbands,nkpts))
  do i=1,nkpts
    read(12,*);read(12,*)
    do j=1,nbands
      read(12,*) temp,eigen(j,i)
    end do
  end do
  close(12)
  allocate(kpath(nkpts))
  call get_kpath(kpath,nkpts,nk)
end subroutine get_band

subroutine get_kpath(kpath,nkpts,nk)
  implicit none
  real(kind=8),parameter :: pi=3.141592653589793D0 
  real(kind=8) :: kpath(nkpts)
  integer :: nkpts,nk,ng
  real(kind=8) :: a(3,3),b(3,3)
  real(kind=8) :: scal,omega
  real(kind=8) :: kpts(3,2),dk,kinit
  integer :: i,j,k

  open(11,file='POSCAR')
  read(11,*)
  read(11,*) scal
  read(11,*) a
  close(11)
  a=a*scal
  omega=dot_product(cross(a(:,1),a(:,2)),a(:,3))
  b(:,1)=2*pi*cross(a(:,2),a(:,3))/omega
  b(:,2)=2*pi*cross(a(:,3),a(:,1))/omega
  b(:,3)=2*pi*cross(a(:,1),a(:,2))/omega
  open(12,file='KPOINTS')
  read(12,*)
  read(12,*) nk
  ng=nkpts/nk
  read(12,*)
  kinit=0.D0
  do i=1,ng
    read(12,*)
    read(12,*) kpts(:,1)
    read(12,*) kpts(:,2)
    kpts=matmul(b,kpts)
    dk=dsqrt(sum((kpts(:,2)-kpts(:,1))**2))/(nk-1)
    forall (j=1:nk)
      kpath((i-1)*nk+j)=kinit+dk*(j-1)
    end forall
    kinit=kpath(i*nk)
  end do
  close(12)
  return
end subroutine get_kpath

function cross(A,B)
  real(kind=8) :: cross(3)
  real(kind=8) :: A(3),B(3)
  cross(1)=A(2)*B(3)-A(3)*B(2)
  cross(2)=A(3)*B(1)-A(1)*B(3)
  cross(3)=A(1)*B(2)-A(2)*B(1)
  return
end function cross

end module eig2band

program effmass
  use eig2band
  implicit none
  real(kind=8),parameter :: pi=3.141592653589793D0
  real(kind=8),parameter :: hbar=6.582119144686551D-16
  real(kind=8),parameter :: me=9.10953D-31
  real(kind=8),parameter :: e=1.602176565D-19
  real(kind=8),allocatable :: dat(:,:)
  real(kind=8),allocatable :: k(:),band(:,:)
  real(kind=8),allocatable :: dk(:),mass(:,:)
  real(kind=8),allocatable :: db(:,:),ddb(:,:)
  real(kind=8),allocatable :: eigen(:,:)
  integer :: nkpts,nk,ng,nb,i

  write(*,'(A,$)') 'Which band do you want to calculate: '
  read(*,*) nb
  write(*,*)
  call get_band(eigen,k,nkpts,nk) 
  ng=nkpts/nk
  allocate(band(nkpts,2))
  allocate(dk(ng))
  allocate(mass(nkpts,2))
  allocate(db(nkpts,2))
  allocate(ddb(nkpts,2))
  band(:,1)=eigen(nb,:)
  band(:,2)=eigen(nb+1,:)
  forall (i=1:ng)
    dk(i)=k(nk*(i-1)+2)-k(nk*(i-1)+1)
  end forall
  do i=1,ng
    db(nk*(i-1)+1:nk*i,1)=gradient(band(nk*(i-1)+1:nk*i,1),dk(i))
    db(nk*(i-1)+1:nk*i,2)=gradient(band(nk*(i-1)+1:nk*i,2),dk(i))
    db(nk*(i-1)+1,:)=0.D0
    db(nk*i,:)=0.D0
    ddb(nk*(i-1)+1:nk*i,1)=gradient(db(nk*(i-1)+1:nk*i,1),dk(i))
    ddb(nk*(i-1)+1:nk*i,2)=gradient(db(nk*(i-1)+1:nk*i,2),dk(i))
  end do
  mass=(1D20*e*hbar**2)/(me*ddb)
  open(11,file='mass.dat')
  write(11,'(A)') '+--------------------------------------------------------------------------+'
  write(11,'(A)') '          K-path     energy(VB)   mass(VB)    energy(CB)   mass(CB)         '
  write(11,'(A)') '          (A^-1)        (eV)        (me)         (eV)        (me)           '
  write(11,'(A)') '+--------------------------------------------------------------------------+'
  do i=1,nkpts
    write(11,'(I6,F11.5,F12.5,Es13.2,F12.5,Es13.2)') i,k(i),band(i,1),mass(i,1),band(i,2),mass(i,2)
    if ( mod(i,nk) == 0 ) then
      write(11,'(A)') '+--------------------------------------------------------------------------+'
    end if
  end do
  close(11)
contains

function gradient(A,d)
  implicit none
  real(kind=8),allocatable :: gradient(:)
  real(kind=8) :: A(:),d
  integer :: n,i

  n=size(A)
  allocate(gradient(n))
  gradient(1)=(A(2)-A(1))/d
  gradient(n)=(A(n)-A(n-1))/d
  forall (i=2:n-1)
    gradient(i)=(A(i+1)-A(i-1))/(2*d)
  end forall
  return
end function gradient

end program effmass
