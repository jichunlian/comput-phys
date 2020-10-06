program hseband
  implicit none
  integer(2) :: nkpts,nibzk,nbands,temp
  real(8),allocatable :: eigen(:,:)
  real(8),allocatable :: kpath(:)
  real(8) :: Efermi
  integer(2) :: i,j
  character(len=4) :: fm

  open(11,file='DOSCAR')
  read(11,*);read(11,*);read(11,*);read(11,*);read(11,*)
  read(11,*) temp,temp,temp,Efermi
  close(11)

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
  eigen=eigen-Efermi

  open(13,file='IBZKPT')
  read(13,*)
  read(13,*) nibzk
  close(13)

  allocate(kpath(nkpts-nibzk))
  call get_kpath(kpath,nkpts-nibzk)

  write(fm,*) nbands
  open(14,file='band.dat')
  do i=nibzk+1,nkpts
    write(14,'(F8.4,'//fm//'F12.6)') kpath(i-nibzk),eigen(:,i)
  end do
  close(14)
end program hseband

subroutine get_kpath(kpath,nkpts)
  implicit none
  real(8),parameter :: pi=3.141592653589793D0 
  real(8) :: kpath(nkpts)
  integer(2) :: nkpts,nk,ng
  real(8) :: a(3,3),b(3,3)
  real(8) :: scal,omega
  real(8) :: kpts(3,2),dk,kinit
  integer(2) :: i,j,k

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
  open(12,file='inkpt')
  read(12,*);read(12,*);read(12,*);read(12,*)
  read(12,*) nk,ng
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
contains

function cross(A,B)
  real(8) :: cross(3)
  real(8) :: A(3),B(3)
  cross(1)=A(2)*B(3)-A(3)*B(2)
  cross(2)=A(3)*B(1)-A(1)*B(3)
  cross(3)=A(1)*B(2)-A(2)*B(1)
  return
end function cross

end subroutine get_kpath
