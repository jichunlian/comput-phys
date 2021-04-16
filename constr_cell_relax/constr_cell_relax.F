!-----------------------------------------------------------------------
!
! At present, VASP DOes not allow to relax the cellshape selectively
! i.e. for instance only cell relaxation in x direction.
! To be more precisse, this behaviour can not be achived via the INCAR
! or POSCAR file.
! However, it is possible to set selected components of the stress tensor
! to zero.
! The most conveninent position to DO this is the routines 
! CONSTR_CELL_RELAX  (constraint cell relaxation).
! FCELL contains the forces on the basis vectors.
! These forces are used to modify the basis vectors according
! to the following equations:
!
!      A_OLD(1:3,1:3)=A(1:3,1:3) ! F90 style 
!      DO J=1,3
!      DO I=1,3
!      DO K=1,3
!        A(I,J)=A(I,J) + FCELL(I,K)*A_OLD(K,J)*STEP_SIZE
!      ENDDO
!      ENDDO
!      ENDDO
! where A holds the basis vectors (in cartesian coordinates).
!
!-----------------------------------------------------------------------
!
!     SUBROUTINE CONSTR_CELL_RELAX(FCELL)
!     USE prec
!     REAL(q) FCELL(3,3)
!
!     just one simple example
!     relaxation in x directions only
!      SAVE=FCELL(1,1)
!      FCELL=0   ! F90 style: set the whole array to zero
!      FCELL(1,1)=SAVE
!     relaxation in z direction only
!      SAVE=FCELL(3,3)
!      FCELL=0   ! F90 style: set the whole array to zero
!      FCELL(3,3)=SAVE
!
!     RETURN
!     END SUBROUTINE
!
!-----------------------------------------------------------------------


! The following code is written by Ji-Chun Lian

SUBROUTINE CONSTR_CELL_RELAX(FCELL,A)
  USE prec
  REAL(q) :: FCELL(3,3),A(3,3)
  LOGICAL(1) :: ALIVE
  INTEGER(1) :: C(9)
  REAL(q) :: F(9),B(9),M(9,9),A0(3,3)

  INQUIRE(FILE='OPTCELL',EXIST=ALIVE)
  IF ( ALIVE == .FALSE. ) RETURN
  OPEN(11,FILE='OPTCELL')
  READ(11,*) C(1),C(4),C(7)
  READ(11,*) C(2),C(5),C(8)
  READ(11,*) C(3),C(6),C(9)
  CLOSE(11)
  F(1:3)=FCELL(1,:)
  F(4:6)=FCELL(2,:)
  F(7:9)=FCELL(3,:)
  M=0.D0
  A0=TRANSPOSE(A)
  M(1:3,1:3)=A0
  M(4:6,4:6)=A0
  M(7:9,7:9)=A0
  B=MATMUL(M,F)
  WHERE ( C == 0 )
    B=0.D0
  END WHERE
  F=MATMUL(INVERSE(M,9),B)
  FCELL(1,:)=F(1:3)
  FCELL(2,:)=F(4:6)
  FCELL(3,:)=F(7:9)
  RETURN

CONTAINS

FUNCTION INVERSE(A,N)
  INTEGER(1) :: N
  REAL(q) :: INVERSE(N,N)
  REAL(q) :: A(N,N)
  REAL(q) :: B(N,N),IA(N,N)
  INTEGER(1) :: I,J

  FORALL(I=1:N,J=1:N,I==J) IA(I,J)=1.D0
  FORALL(I=1:N,J=1:N,I/=J) IA(I,J)=0.D0
  B=A
  CALL UPPER(B,IA,N)
  CALL LOWER(B,IA,N)
  FORALL(I=1:N) IA(I,:)=IA(I,:)/B(I,I)
  INVERSE=IA
  RETURN
END FUNCTION INVERSE

SUBROUTINE UPPER(M,S,N)
  INTEGER(1) :: N
  REAL(q) :: M(N,N)
  REAL(q) :: S(N,N)
  INTEGER(1) :: I,J
  REAL(q) :: E

  DO I=1,N-1
    DO J=I+1,N
      E=M(J,I)/M(I,I)
      M(J,I:N)=M(J,I:N)-M(I,I:N)*E
      S(J,:)=S(J,:)-S(I,:)*E
    END DO
  END DO
  RETURN
END SUBROUTINE UPPER

SUBROUTINE LOWER(M,S,N)
  INTEGER(1) :: N
  REAL(q) :: M(N,N)
  REAL(q) :: S(N,N)
  INTEGER(1) :: I,J
  REAL(q) :: E

  DO I=N,2,-1
    DO J=I-1,1,-1
      E=M(J,I)/M(I,I)
      M(J,I:N)=M(J,I:N)-M(I,I:N)*E
      S(J,:)=S(J,:)-S(I,:)*E
    END DO
  END DO
  RETURN
END SUBROUTINE LOWER

END SUBROUTINE CONSTR_CELL_RELAX
