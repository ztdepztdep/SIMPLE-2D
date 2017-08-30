  
  SUBROUTINE PRESSURE
  USE GLOBAL
  IMPLICIT NONE   
  INTEGER(KIND=4)::INCOUNTER
  REAL(KIND=4)::R0,RK,P1,P2
    NF=3	  
	CON=0	  
IF(LSOLVE(NF) .eqv. .TRUE.) THEN
    IST=2
    JST=2
    ISTF=IST-1
    JSTF=JST-1
    CSUM=0.0
	CALL GAMSOR

DO I=2,L2
	ARHO=R(1)*XCV(I)*RHO(I,1)    
	CON(I,2)=CON(I,2)+ARHO*MIV(I,2)
    AJM(I,2)=0.0
ENDDO   

DO J=2,M2
    ARHO=ARX(J)*RHO(1,J)
    CON(2,J)=CON(2,J)+ARHO*MIU(2,J)
    AIM(2,J)=0. 
DO I=2,L2
  IF(I .EQ. L2) THEN
	ARHO=ARX(J)*RHO(L1,J)
    CON(I,J)=CON(I,J)-ARHO*MIU(L1,J)
    AIP(I,J)=0.0
  ELSE
    ARHO=ARX(J)*( FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J) )
    FLOW=ARHO*MIU(I+1,J)
    CON(I,J)=CON(I,J)-FLOW
    CON(I+1,J)=CON(I+1,J)+FLOW
    AIP(I,J)=ARHO*DEU(I+1,J)
    AIM(I+1,J)=AIP(I,J)
  ENDIF
  IF(J .EQ. M2) THEN
    ARHO=RMN(M1)*XCV(I)*RHO(I,M1)
    CON(I,J)=CON(I,J)-ARHO*MIV(I,M1)
   	AJP(I,J)=0.
  ELSE   
    ARHO=RMN(J+1)*XCV(I)*( FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J) )
    FLOW=ARHO*MIV(I,J+1)
    CON(I,J)=CON(I,J)-FLOW
    CON(I,J+1)=CON(I,J+1)+FLOW
	
    AJP(I,J)=ARHO*DNV(I,J+1)
    AJM(I,J+1)=AJP(I,J)
  ENDIF
    AP(I,J)=AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J)      
	PC(I,J)=0.0	      
	CSUM=MAX( CSUM,ABS( CON(I,J) ) )
  ENDDO
  ENDDO
    CSUM=CSUM/REFDENOM

    CALL SOLVER

DO J=2,M2
DO I=2,L2
    IF(I .NE. 2) THEN
	MIU(I,J)=MIU(I,J)+DEU(I,J)*( PC(I-1,J)-PC(I,J) ) 
	ENDIF

    IF(J .NE. 2) THEN
	MIV(I,J)=MIV(I,J)+DNV(I,J)*( PC(I,J-1)-PC(I,J) )
	ENDIF
ENDDO
ENDDO
  	  P=P+0.3*PC

DO J=2,M2
DO I=2,L2
	P1= FX(I)*PC(I,J)+FXM(I)*PC(I-1,J)  
    P2= FX(I+1)*PC(I+1,J)+FXM(I+1)*PC(I,J)
    U(I,J)=U(I,J)+DU(I,J)*( P1-P2 ) 
	
	P1= FY(J)*PC(I,J)+FYM(J)*PC(I,J-1)
    P2= FY(J+1)*PC(I,J+1)+FYM(J+1)*PC(I,J)
    V(I,J)=V(I,J)+DV(I,J)*( P1-P2 )
ENDDO
ENDDO

ENDIF
      
	  STREAM(2,2)=0
   	  DO I=2,L1
	  IF(I.NE.2)THEN
	  STREAM(I,2)=STREAM(I-1,2)-V(I-1,2)*R(1)*XCV(I-1)
	  ENDIF
	  DO J=3,M1
	  STREAM(I,J)=STREAM(I,J-1)+U(I,J-1)*ARX(J-1)
	  ENDDO
	  ENDDO

      END


  SUBROUTINE UPDATEPRESSURE
  USE GLOBAL
  IMPLICIT NONE   
	  DO J=2,M2
       P(1,J)=P(2,J)+XDIF(2)*( P(2,J)-P(3,J) )/XDIF(3)
       P(L1,J)=P(L2,J)+XDIF(L1)*( P(L2,J)-P(L3,J) )/XDIF(L2)
      ENDDO
      
	  DO I=2,L2
       P(I,1)=P(I,2)+YDIF(2)*( P(I,2)-P(I,3) )/YDIF(3)
       P(I,M1)=P(I,M2)+YDIF(M1)*( P(I,M2)-P(I,M3) )/YDIF(M2)
      ENDDO

      P(1,1)=P(2,1)+P(1,2)-P(2,2)
      P(L1,1)=P(L2,1)+P(L1,2)-P(L2,2)
      P(1,M1)=P(2,M1)+P(1,M2)-P(2,M2)
      P(L1,M1)=P(L2,M1)+P(L1,M2)-P(L2,M2)
      P=P-P(1,1)
END
