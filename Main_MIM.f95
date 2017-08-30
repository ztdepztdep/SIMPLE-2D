  SUBROUTINE MIM
  USE GLOBAL
  IMPLICIT NONE   
REAL*8 P1,P2,TMEP

 DO J=2,M2
 DO I=2,L2
   UHAT(I,J)=( UAIP(I,J)*U(I+1,J)+UAIM(I,J)*U(I-1,J)+UAJP(I,J)*U(I,J+1)+UAJM(I,J)*U(I,J-1)+UCON(I,J) )/UAP(I,J)
   VHAT(I,J)=( VAIP(I,J)*V(I+1,J)+VAIM(I,J)*V(I-1,J)+VAJP(I,J)*V(I,J+1)+VAJM(I,J)*V(I,J-1)+VCON(I,J) )/VAP(I,J)
 ENDDO
 ENDDO 



 DO J=2,M2
 DO I=3,L2
   DEU(I,J)=FX(I)*DU(I,J)+FXM(I)*DU(I-1,J)
   TEMP=FX(I)*UHAT(I,J)+FXM(I)*UHAT(I-1,J)+DEU(I,J)*( P(I-1,J)-P(I,J) )
   MIU(I,J)=TEMP+1.0/(1.+EU)*MIU(I,J)

   !TEMP=FX(I)*UAP(I,J)+FXM(I)*UAP(I-1,J)
   !VOL=YCVR(J)*XDIF(I) 
   !P1=( P(I,J)-P(I-1,J) )/XDIF(I)
   !MIU(I,J)=FX(I)*U(I,J)+FXM(I)*U(I-1,J)-VOL*TEMP*( P1-0.5*( DPX(I,J)+DPX(I-1,J) ) )
 ENDDO
 ENDDO 



 DO J=3,M2
 DO I=2,L2
   DNV(I,J)=FY(J)*DV(I,J)+FYM(J)*DV(I,J-1)
   TEMP=FY(J)*VHAT(I,J)+FYM(J)*VHAT(I,J-1)+DNV(I,J)*( P(I,J-1)-P(I,J) )
   MIV(I,J)=TEMP+1.0/(1.+EU)*MIV(I,J)


  ! TEMP=FY(J)*VAP(I,J)+FYM(J)*VAP(I,J-1)
  ! VOL=YCVR(J)*XDIF(I) 
  ! P1=( P(I,J)-P(I,J-1) )/YDIF(J)
   !MIV(I,J)=FY(J)*V(I,J)+FYM(J)*V(I,J-1)-VOL*TEMP*( P1-0.5*( DPY(I,J)+DPY(I,J-1) ) )
 ENDDO
 ENDDO 

 
 END