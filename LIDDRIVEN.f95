

SUBROUTINE DEFINE
USE GLOBAL 
IMPLICIT NONE 
   
  Mode=1

  LSOLVE(1)=.false.
  LSOLVE(2)=.false.
  LSOLVE(3)=.false.
  LSOLVE(5)=.true.

       
  NTIMES(1)=1
  NTIMES(2)=1
  NTIMES(3)=3
  NTIMES(5)=2
  DT=1E+30

  RHOCON=1.0
  TOL=1E-8

  RENOLD=1000.0

 EU=4
 EOTHER=5.0
End

SUBROUTINE UGRID
USE GLOBAL
IMPLICIT NONE
 XL=1.0
 YL=1.0
END

SUBROUTINE START
USE GLOBAL
IMPLICIT NONE
U=0
V=0
U(:,M1)=0.0
T=0.0

T(1,:)=1.0
T(L1,:)=0.0

  MIU(2,:)=U(1,:)
  MIU(L1,:)=U(L1,:)
  MIU(:,1)=U(:,1)
  MIU(:,M1)=U(:,M1)
 
  MIV(:,2)=V(:,1)
  MIV(:,M1)=V(:,M1)
  MIV(1,:)=V(1,:)
  MIV(L1,:)=V(L1,:)
END

SUBROUTINE	 BOUND
USE GLOBAL 
IMPLICIT NONE
REAL(KIND=4)::F
MONITOR(1)=U(1,1)
MONITOR(2)=T(3,3)

T(:,M1)=T(:,M2)
T(:,1)=T(:,2)
END	

SUBROUTINE GAMSOR
USE GLOBAL 
IMPLICIT NONE

IF(NF==1) THEN
sp=0

SC=0
GAM=1.0/RENOLD
ENDIF

IF(NF==2) THEN
sp=0

SC=0
GAM=1.0/RENOLD
ENDIF

IF(NF==5) THEN
gam=1.0
gam(:,M1)=0.0
gam(:,1)=0.0

ENDIF
END

SUBROUTINE DENSE
USE GLOBAL
IMPLICIT NONE
RHO=RHOCON
END
 

SUBROUTINE OUTPUT	!FOR TECPLOT 
USE GLOBAL 
IMPLICIT NONE

	OUTNAME(1)="STREAMLINE.dat"
    OPEN(UNIT=10,FILE=OUTNAME(1))
    WRITE(10,*) 'VARIABLES="X","Y","STREAM LINE"'
    WRITE(10,*) 'ZONE I=',L2,',J=',M2,'F=POINT'
     DO J=2,M1
     DO I=2,L1
   	 WRITE(10,*)  XU(I),YV(J),STREAM(I,J)
     ENDDO
    ENDDO
	CLOSE(10)

	WRITE(FILENAME,"(f7.3)") RENOLD
	OUTNAME(1)="UV.dat"
    OPEN(UNIT=10,FILE=OUTNAME(1))
    WRITE(10,*) 'VARIABLES="X","Y","U","V","T"'
    WRITE(10,*) 'ZONE I=',L1,',J=',M1,'F=POINT'
     DO J=1,M1
     DO I=1,L1
   	 WRITE(10,*)  x(I),y(J),U(I,J),V(I,J),T(I,J)
     ENDDO
    ENDDO
	CLOSE(10)
   
   OPEN(UNIT=10,FILE='centeru.dat')
     DO J=1,M1
      WRITE(10,'(3f10.4)') y(J),MIU(L1/2+1,J),0.5*( U( L1/2+1,J )+U( L1/2,J) )
     ENDDO
    close(10) 

	   OPEN(UNIT=10,FILE='centerv.dat')
     DO I=1,L1
      WRITE(10,'(3f10.4)') x(i),MIV(I,M1/2+1),0.5*( V(I,M1/2+1)+V(I,M1/2)  )
     ENDDO
      close(10)

END

