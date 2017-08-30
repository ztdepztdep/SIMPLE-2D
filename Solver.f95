SUBROUTINE QUICK0
  USE GLOBAL
  IMPLICIT NONE

 IF( I>IST.AND.I<L2.AND.J>JST.AND.J<M2) THEN

CON(I,J)=CON(I,J)+MAX(FLUXE(I,J),0.0 )*(  PHI(I-1,J,NF)/8.0+2*PHI(I,J,NF)/8.0-3*PHI(I+1,J,NF)/8.0  )&
                &-MAX(-FLUXE(I,J),0.0)*(  -3*PHI(I,J,NF)/8.0+2*PHI(I+1,J,NF)/8.0+PHI(I+2,J,NF)/8.0 )
				!SADE
CON(I,J)=CON(I,J)-(  MAX(FLUXE(I-1,J),0.0 )*(  PHI(I-2,J,NF)/8.0+2*PHI(I-1,J,NF)/8.0-3*PHI(I,J,NF)/8.0  )&
                &-MAX(-FLUXE(I-1,J),0.0)*(  -3*PHI(I-1,J,NF)/8.0+2*PHI(I,J,NF)/8.0+PHI(I+1,J,NF)/8.0 )  )
			   	!SADW
CON(I,J)=CON(I,J)+MAX(FLUXN(I,J),0.0 )*(  PHI(I,J-1,NF)/8.0+2*PHI(I,J,NF)/8.0-3*PHI(I,J+1,NF)/8.0  )&
                &-MAX(-FLUXN(I,J),0.0)*(  -3*PHI(I,J,NF)/8.0+2*PHI(I,J+1,NF)/8.0+PHI(I,J+2,NF)/8.0 )
				!SADN
CON(I,J)=CON(I,J)-(  MAX(FLUXN(I,J-1),0.0 )*(  PHI(I,J-2,NF)/8.0+2*PHI(I,J-1,NF)/8.0-3*PHI(I,J,NF)/8.0  )&
                &-MAX(-FLUXN(I,J-1),0.0)*(  -3*PHI(I,J-1,NF)/8.0+2*PHI(I,J,NF)/8.0+PHI(I,J+1,NF)/8.0 )  )
			   !SADS
   ENDIF

if(i==Ist) CON(I,J)=CON(I,J)+MAX(FLUXE(I,J),0.0 )*(  PHI(I-1,J,NF)/8.0+2*PHI(I,J,NF)/8.0-3*PHI(I+1,J,NF)/8.0  )&
                &-MAX(-FLUXE(I,J),0.0)*(  -3*PHI(I,J,NF)/8.0+2*PHI(I+1,J,NF)/8.0+PHI(I+2,J,NF)/8.0 )
				!SADE
if(i==L2)CON(I,J)=CON(I,J)-(  MAX(FLUXE(I-1,J),0.0 )*(  PHI(I-2,J,NF)/8.0+2*PHI(I-1,J,NF)/8.0-3*PHI(I,J,NF)/8.0  )&
                &-MAX(-FLUXE(I-1,J),0.0)*(  -3*PHI(I-1,J,NF)/8.0+2*PHI(I,J,NF)/8.0+PHI(I+1,J,NF)/8.0 )  )
			   	!SADW
if(J==Jst)CON(I,J)=CON(I,J)+MAX(FLUXN(I,J),0.0 )*(  PHI(I,J-1,NF)/8.0+2*PHI(I,J,NF)/8.0-3*PHI(I,J+1,NF)/8.0  )&
                &-MAX(-FLUXN(I,J),0.0)*(  -3*PHI(I,J,NF)/8.0+2*PHI(I,J+1,NF)/8.0+PHI(I,J+2,NF)/8.0 )
				!SADN
if(J==M2)CON(I,J)=CON(I,J)-(  MAX(FLUXN(I,J-1),0.0 )*(  PHI(I,J-2,NF)/8.0+2*PHI(I,J-1,NF)/8.0-3*PHI(I,J,NF)/8.0  )&
                &-MAX(-FLUXN(I,J-1),0.0)*(  -3*PHI(I,J-1,NF)/8.0+2*PHI(I,J,NF)/8.0+PHI(I,J+1,NF)/8.0 )  )
			   !SADS
END

SUBROUTINE SOLVER0
USE GLOBAL
IMPLICIT NONE
ISTF=IST-1
JSTF=JST-1

IT1=L2+IST
IT2=L3+IST

JT1=M2+JST
JT2=M3+JST

DO K=1,NTIMES(NF)

DO J=JST,M2
 PT(ISTF)=0.0
 QT(ISTF)=PHI(ISTF,J,NF) 
 DO I=IST,L2
   DENOM=AP(I,J)-PT(I-1)*AIM(I,J)+1E-30
   PT(I)=AIP(I,J)/DENOM
   TEMP=CON(I,J)+AJP(I,J)*PHI(I,J+1,NF)+AJM(I,J)*PHI(I,J-1,NF)
   QT(I)=( TEMP+AIM(I,J )*QT(I-1) )/DENOM
 ENDDO 
 DO II=IST,L2
  I=IT1-II
  PHI(I,J,NF)=PHI(I+1,J,NF)*PT(I)+QT(I)
 ENDDO
ENDDO

DO JJ=JST,M3
J=JT2-JJ
 PT(ISTF)=0.0
 QT(ISTF)=PHI(ISTF,J,NF)
 DO I=IST,L2
   DENOM=AP(I,J)-PT(I-1)*AIM(I,J)+1E-30
   PT(I)=AIP(I,J)/DENOM
   TEMP=CON(I,J)+AJP(I,J)*PHI(I,J+1,NF)+AJM(I,J)*PHI(I,J-1,NF)
   QT(I)=(TEMP+AIM(I,J)*QT(I-1))/DENOM
 ENDDO 
 DO II=IST,L2
  I=IT1-II
  PHI(I,J,NF)=PHI(I+1,J,NF)*PT(I)+QT(I)
 ENDDO
ENDDO

DO I=IST,L2
 PT(JSTF)=0.0
 QT(JSTF)=PHI(I,JSTF,NF)
 DO J=JST,M2
  DENOM=AP(I,J)-PT(J-1)*AJM(I,J)+1E-30
  PT(J)=AJP(I,J)/(DENOM)
  TEMP=CON(I,J)+AIP(I,J)*PHI(I+1,J,NF)+AIM(I,J)*PHI(I-1,J,NF)
  QT(J)=(TEMP+AJM(I,J)*QT(J-1))/DENOM
 ENDDO
 DO JJ=JST,M2
 J=JT1-JJ
 PHI(I,J,NF)=PHI(I,J+1,NF)*PT(J)+QT(J)
  ENDDO
ENDDO

DO II=IST,L3
 I=IT2-II
 PT(JSTF)=0.0
 QT(JSTF)=PHI(I,JSTF,NF)
 DO J=JST,M2
  DENOM=AP(I,J)-PT(J-1)*AJM(I,J)+1E-30
  PT(J)=AJP(I,J)/DENOM
  TEMP=CON(I,J)+AIP(I,J)*PHI(I+1,J,NF)+AIM(I,J)*PHI(I-1,J,NF)
  QT(J)=(TEMP+AJM(I,J)*QT(J-1))/DENOM
 ENDDO
  DO JJ=JST,M2
  J=JT1-JJ
  PHI(I,J,NF)=PHI(I,J+1,NF)*PT(J)+QT(J)
  ENDDO
ENDDO
ENDDO
END




SUBROUTINE SOLVER
USE GLOBAL
IMPLICIT NONE
REAL(KIND=8)::LW(L1,M1),LS(L1,M1),LPR(L1,M1),UUN(L1,M1), UUE(L1,M1),RES(L1,M1)

REAL(KIND=8)::P1,P2,RESL,ALFA=0.9,RESNOR(NMAX),RSM
ISTF=IST-1
JSTF=JST-1
IT1=L2+IST
IT2=L3+IST
JT1=M2+JST
JT2=M3+JST
      DO I=IST,L2
        DO J=JST,M2
          LW(I,J)=-AIM(I,J)/(1.+ALFA*UUN(I-1,J))
          LS(I,J)=-AJM(I,J)/(1.+ALFA*UUE(I,J-1))
          P1=ALFA*LW(I,J)*UUN(I-1,J)
          P2=ALFA*LS(I,J)*UUE(I,J-1)
          LPR(I,J)=1./(AP(I,J)+P1+P2-LW(I,J)*UUE(I-1,J)-LS(I,J)*UUN(I,J-1)+1.E-30)
          UUN(I,J)=(-AJP(I,J)-P1)*LPR(I,J)
          UUE(I,J)=(-AIP(I,J)-P2)*LPR(I,J)
        END DO
      END DO

      DO K=1,NTIMES(NF)
	    RESL=0.
        DO I=IST,L2
          DO J=JST,M2
            RES(I,J)=CON(I,J)-AP(I,J)*PHI(I,J,NF)+AJP(I,J)*PHI(I,J+1,NF)+&
                     AJM(I,J)*PHI(I,J-1,NF)+AIP(I,J)*PHI(I+1,J,NF)+AIM(I,J)*PHI(I-1,J,NF)
            RESL=RESL+ABS( RES(I,J) )
            RES(I,J)=( RES(I,J)-LS(I,J)*RES(I,J-1)-LW(I,J)*RES(I-1,J))*LPR(I,J)
          END DO
        END DO

        IF(K.EQ.1) RESNOR(NF)=RESL
        RSM=RESL/RESNOR(NF)

        DO I=L2,IST,-1
          DO J=M2,JST,-1
            RES(I,J)=RES(I,J)-UUN(I,J)*RES(I,J+1)-UUE(I,J)*RES(I+1,J)
            PHI(I,J,NF)=PHI(I,J,NF)+RES(I,J)
          END DO
        END DO
!.....CHECK CONVERGENCE
        
     !  IF(RSM<0.2) RETURN

      END DO
      END




	 SUBROUTINE SOLVER3 ! CGSTAB
      USE GLOBAL
	  IMPLICIT NONE
     REAL*8::RESMAX,ALF,ERRNOR,RES0,GAAM,BET,BETO,OM,UKRESO,RESL,VRES,VV,RSM
     REAL*8::PK(L1,M1),ZK(L1,M1),D(L1,M1),RES(L1,M1),RESO(L1,M1),UK(L1,M1),VK(L1,M1)
	 ISTF=IST-1
	 JSTF=JST-1
	 IT1=L2+IST
	 IT2=L3+IST
	 JT1=M2+JST
	 JT2=M3+JST
!.....CALCULATE INITIAL RESIDUAL VECTOR
          RES0=0.0
          DO I=IST,L2
          DO J=JST,M2 
          RES(I,J)=CON(I,J)-AP(I,J)*PHI(I,J,NF)+AIP(I,J)*PHI(I+1,J,NF)+AIM(I,J)*PHI(I-1,J,NF)+&
                &AJP(I,J)*PHI(I,J+1,NF)+AJM(I,J)*PHI(I,J-1,NF)
          RES0=RES0+ABS( RES(I,J) )
         END DO
         END DO	     
!.....CALCULATE ELEMENTS OF PRECONDITIONING MATRIX DIAGONAL
          DO I=IST,L2
          DO J=JST,M2
          D(I,J)=1.0/(1E-30+AP(I,J)-AIM(I,J)*D(I-1,J)*AIP(I-1,J)-AJM(I,J)*D(I,J-1)*AJP(I,J-1) ) 
          END DO
          END DO
!.....INITIALIZE WORKING ARRAYS AND CONSTANTS
          DO I=IST,L2
          DO J=JST,M2
            RESO(I,J)=RES(I,J)
            PK(I,J)=0.
            UK(I,J)=0.
            ZK(I,J)=0.
            VK(I,J)=0.
          END DO
        END DO
        ALF=1.
        BETO=1.
        GAAM=1.

!.....START INNER ITERATIONS

      DO K=1,NTIMES(NF)
!..... CALCULATE BETA AND OMEGA
          BET=0.
          DO I=IST,L2
          DO J=JST,M2
           BET=BET+RES(I,J)*RESO(I,J)
          END DO
          END DO
           OM=BET*GAAM/( ALF*BETO+1.E-30 )
           BETO=BET
!..... CALCULATE PK
        DO I=IST,L2
          DO J=JST,M2
            PK(I,J)=RES(I,J)+OM*( PK(I,J)-ALF*UK(I,J) )
          END DO
        END DO
!.....SOLVE (M ZK = PK) - FORWARD SUBSTITUTION
        DO I=IST,L2
          DO J=JST,M2            
            ZK(I,J)=( PK(I,J)+AIM(I,J)*ZK(I-1,J)+AJM(I,J)*ZK(I,J-1) )*D(I,J)
          END DO
        END DO

          DO I=IST,L2
          DO J=JST,M2
            ZK(I,J)=ZK(I,J)/( D(I,J)+1.E-30 )
          END DO
        END DO
!..... BACKWARD SUBSTITUTION

        DO I=L2,IST,-1
          DO J=M2,JST,-1            
            ZK(I,J)=( ZK(I,J)+AIP(I,J)*ZK(I+1,J)+AJP(I,J)*ZK(I,J+1) )*D(I,J)
          END DO
        END DO
!.....CALCULATE UK = A.PK
        DO I=IST,L2
          DO J=JST,M2            
            UK(I,J)=AP(I,J)*ZK(I,J)-AIP(I,J)*ZK(I+1,J)-&
                    &AIM(I,J)*ZK(I-1,J)-AJP(I,J)*ZK(I,J+1)-&
                    &AJM(I,J)*ZK(I,J-1)
          END DO
        END DO
!..... CALCULATE SCALAR PRODUCT UK.RESO AND GAMMA
       UKRESO=0.
        DO I=IST,L2
          DO J=JST,M2
            UKRESO=UKRESO+UK(I,J)*RESO(I,J)
          END DO
        END DO
      GAAM=BET/(1E-30+UKRESO)
!.....UPDATE (FI) AND CALCULATE W (OVERWRITE RES - IT IS RES-UPDATE)
        DO I=IST,L2
          DO J=JST,M2
            PHI(I,J,NF)=PHI(I,J,NF)+GAAM*ZK(I,J)
            RES(I,J)=RES(I,J)-GAAM*UK(I,J)
          END DO
        END DO
!.....SOLVE (M Y = W); Y OVERWRITES ZK; FORWARD SUBSTITUTION
        DO I=IST,L2
          DO J=JST,M2
            ZK(I,J)=( RES(I,J)+AIM(I,J)*ZK(I-1,J)+AJM(I,J)*ZK(I,J-1) )*D(I,J)
           END DO
         END DO

        DO I=IST,L2
          DO J=JST,M2          
            ZK(I,J)=ZK(I,J)/( D(I,J)+1.E-30 )
          END DO
        END DO

!.....BACKWARD SUBSTITUTION

        DO I=L2,IST,-1
          DO J=M2,JST,-1
            ZK(I,J)=( ZK(I,J)+AIP(I,J)*ZK(I+1,J)+AJP(I,J)*ZK(I,J+1) )*D(I,J)
          END DO
        END DO
!.....CALCULATE V = A.Y (VK = A.ZK)
        DO I=IST,L2
          DO J=JST,M2
            VK(I,J)=AP(I,J)*ZK(I,J)-AIP(I,J)*ZK(I+1,J)-&
                   &AIM(I,J)*ZK(I-1,J)-AJP(I,J)*ZK(I,J+1)-&
                   &AJM(I,J)*ZK(I,J-1) 
          END DO
        END DO
!..... CALCULATE ALPHA (ALF)
      VRES=0.
      VV=0.
        DO I=IST,L2
          DO J=JST,M2
            VRES=VRES+VK(I,J)*RES(I,J)
            VV=VV+VK(I,J)*VK(I,J)
          END DO
      END DO

      ALF=VRES/(VV+1.E-30)
!.....UPDATE VARIABLE (FI) AND RESIDUAL (RES) VECTORS
        RESL=0.
        DO I=IST,L2
          DO J=JST,M2
            PHI(I,J,NF)=PHI(I,J,NF)+ALF*ZK(I,J)
            RES(I,J)=RES(I,J)-ALF*VK(I,J)
            RESL=RESL+ABS( RES(I,J) )
          END DO
      END DO

	ENDDO
    END
