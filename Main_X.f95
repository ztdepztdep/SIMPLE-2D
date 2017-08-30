
SUBROUTINE XVELOCITY
USE GLOBAL
IMPLICIT NONE
REAL*8 P1,P2
 NF=1
 IF(LSOLVE(NF).eqv..TRUE.) THEN 
 IST=2
 JST=2 
 ISTF=IST-1
 JSTF=JST-1
 CALL GAMSOR

DO I=2,L2
     AREA=R(1)*XCV(I)    
	 FLOW=AREA*MIV(I,2)*RHO(I,1)     
	 DIF=AREA*GAM(I,1)/YDIF(2) 
	 AJM(I,2)=DIF+MAX(0.0,FLOW)
ENDDO

DO J=2,M2
     FLOW=ARX(J)*MIU(2,J)*RHO(1,J)
     DIF=ARX(J)*GAM(1,J)/(XDIF(2)*SX(J))
     AIM(2,J)=DIF+MAX(0.0,FLOW)  
DO I=2,L2  
  IF(I==L2) THEN 
      FLOW=ARX(J)*MIU(L1,J)*RHO(L1,J)
      DIF=ARX(J)*GAM(L1,J)/( XDIF(L1)*SX(J) )
  ELSE
	  FLOW=ARX(J)*MIU(I+1,J)*( FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J) )
      DIF=ARX(J)*2.0*GAM(I,J)*GAM(I+1,J)/( (XCV(I)*GAM(I+1,J)+XCV(I+1)*GAM(I,J)+1.0E-30)*SX(J) ) 
  ENDIF    
	  AIM(I+1,J)=DIF+MAX(0.,FLOW)
      AIP(I,J)=AIM(I+1,J)-FLOW 	    
	  FLUXE(I,J)=FLOW 

	  AREA=RMN(J+1)*XCV(I)	
   IF(J==M2) THEN
      FLOW=AREA*MIV(I,M1)*RHO(I,M1)
      DIF=AREA*GAM(I,M1)/YDIF(M1)
   ELSE
	  FLOW=AREA*MIV(I,J+1)*( FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J) )
      DIF=AREA*2.*GAM(I,J)*GAM(I,J+1)/( YCV(J)*GAM(I,J+1)+YCV(J+1)*GAM(I,J)+1.0E-30 )	
   ENDIF 
	  AJM(I,J+1)=DIF+MAX(0.,FLOW)    
	  AJP(I,J)=AJM(I,J+1)-FLOW   
	  FLUXN(I,J)=FLOW 


      VOL=YCVR(J)*XCV(I)   
	  APT=RHO(I,J)/DT    
	  AP(I,J)=SP(I,J)-APT
	  CON(I,J)=SC(I,J)+APT*U(I,J)
	  AP(I,J)=( -AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J) )

	  CON(I,J)=CON(I,J)*VOL

      CALL QUICK0
	  UCON(I,J)=CON(I,J)   

   	  AP(I,J)=AP(I,J)*( 1.0+1.0/EU )  
	  CON(I,J)=CON(I,J)+AP(I,J)*U(I,J)/(1+EU)

      DU(I,J)=VOL/( XCV(I)*SX(J) )  
	  DPX(I,J)=(P2-P1)/XCV(I)	 

	  P1= FX(I)*P(I,J)+FXM(I)*P(I-1,J)  
      P2= FX(I+1)*P(I+1,J)+FXM(I+1)*P(I,J)
  
      CON(I,J)=CON(I,J)+DU(I,J)*( P1-P2 ) 
      DU(I,J)=DU(I,J)/ AP(I,J)
 ENDDO
 ENDDO 

   UAP=AP
   UAIM=AIM
   UAIP=AIP
   UAJM=AJM
   UAJP=AJP
  
  
   TEMP=1E-35
   TEMP1=1E-35   
   DO I=IST,L2
   DO J=JST,M2	 
   TEMP=TEMP+(AP(I,J)*PHI(I,J,1)-AIM(I,J)*PHI(I-1,J,1)-AIP(I,J)*PHI(I+1,J,1)-&
              &AJM(I,J)*PHI(I,J-1,1)-AJP(I,J)*PHI(I,J+1,1)-CON(I,J))**2
   TEMP1=TEMP1+(AP(I,J)*PHI(I,J,1))**2
   ENDDO
   ENDDO     
    
   USUM=SQRT(TEMP/TEMP1)
   
  
   I=NX/2
   REFDENOM=1E-30	
   DO J=JST,M2
   REFDENOM=REFDENOM+RHO(I,J)*ABS(PHI(I,J,1))*ARX(J)
   ENDDO
   ENDIF

 CALL SOLVER



   END	

