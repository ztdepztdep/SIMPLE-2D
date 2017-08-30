

SUBROUTINE GRIDSETUP
  USE GLOBAL
  IMPLICIT NONE
  CALL UGRID
  
   XU(2)=0

  DX=XL/REAL(NX)
  DX=1.0/REAL(NX)
  DO I=3,L1
  XU(I)=XU(I-1)+DX
  ENDDO
  
   YV(2)=0
 DY=YL/REAL(NY)
  DY=1.0/REAL(NY) !引入了网格控制函数的重要变化 
  
  DO J=3,M1
  YV(J)=YV(J-1)+DY
  ENDDO

  DO I=2,L1
   XU(I)=XL*( 0.5+0.5*TANH( KEXIX*( 2*XU(I)-1 ) )/TANH( KEXIX ) )

  ENDDO
    XU(L1)=XL
  DO J=2,M1
   YV(J)=YL*( 0.5+0.5*TANH( KEXIY*( 2*YV(J)-1 ) )/TANH( KEXIY ) )
  ENDDO
   YV(M1)=YL

  X(1)=XU(2)
  DO I=2,L2
  X(I)=0.5*(XU(I+1)+XU(I))
  ENDDO
  X(L1)=XU(L1)

  Y(1)=YV(2)
  DO J=2,M2
  Y(J)=0.5*( YV(J+1)+YV(J) )
  ENDDO
  Y(M1)=YV(M1)

  DO I=2,L1
  XDIF(I)=X(I)-X(I-1)
  ENDDO
  XDIF(1)=XDIF(2)

  DO I=2,L2
  XCV(I)=XU(I+1)-XU(I)
  ENDDO
  XCV(1)=0.0
  XCV(L1)=0.0

  DO I=3,L2
  XCVS(I)=XDIF(I)
  ENDDO
  XCVS(3)=XCVS(3)+XDIF(2)
  XCVS(L2)=XCVS(L2)+XDIF(L1)

  DO I=3,L3
  XCVI(I)=0.5*XCV(I)
  XCVIP(I)=XCVI(I)
  ENDDO
  XCVI(L2)=XCV(L2)
  XCVIP(2)=XCV(2)

  DO J=2,M1
  YDIF(J)=Y(J)-Y(J-1)
  ENDDO
  DO J=2,M2
  YCV(J)=YV(J+1)-YV(J)
  ENDDO
  YCV(1)=0.0
  YCV(M1)=0.0

  DO J=3,M2
  YCVS(J)=YDIF(J)
  ENDDO
  YCVS(3)=YCVS(3)+YDIF(2)
  YCVS(M2)=YCVS(M2)+YDIF(M1)
 
  IF(MODE==1) THEN	  
  RMN=1.0
  R=1.0
  ELSE
  DO J=2,M1
  R(J)=R(J-1)+YDIF(J)
  ENDDO
  RMN(2)=R(1)
  DO J=3,M2
  RMN(J)=RMN(J-1)+YCV(J-1)
  ENDDO
  RMN(M1)=R(M1)
  ENDIF

  DO J=1,M1
  SX(J)=1.
  SXMN(J)=1.    
  IF(MODE ==3) THEN
  SX(J)=R(J)
  IF(J .NE. 1) THEN
  SXMN(J)=RMN(J)
  ENDIF
  ENDIF
  ENDDO    

  DO J=2,M2
  YCVR(J)=R(J)*YCV(J)	  	
  ARX(J)=YCVR(J)      
  IF (MODE ==3) THEN
  ARX(J)=YCV(J)
  ENDIF
  ENDDO

  DO J=4,M3
  YCVRS(J)=0.5*( R(J)+R(J-1) )*YDIF(J)
  ENDDO
  YCVRS(3)=0.5*( R(3)+R(1) )*YCVS(3)
  YCVRS(M2)=0.5*( R(M1)+R(M3) )*YCVS(M2)    !计算体积是用，其他用途并无   

  IF(MODE==2) THEN
   DO J=3,M3
   ARXJ(J)=0.25*( 1.+RMN(J)/R(J) )*ARX(J)   
   ARXJP(J)=ARX(J)-ARXJ(J)  
   ENDDO 
   ELSE
   DO J=3,M3
   ARXJ(J)=0.5*ARX(J)
   ARXJP(J)=ARXJ(J)
   ENDDO
  ENDIF
   ARXJP(2)=ARX(2)
   ARXJ(M2)=ARX(M2)
     	 
  DO J=3,M3
  FV(J)=ARXJP(J)/ARX(J)
  FVP(J)=1.-FV(J)          
  ENDDO

  DO I=3,L2
  FX(I)=0.5*XCV(I-1)/XDIF(I)
  FXM(I)=1.-FX(I)
  ENDDO
  FX(2)=0.
  FXM(2)=1.
  FX(L1)=1.
  FXM(L1)=0.

  DO J=3,M2
  FY(J)=0.5*YCV(J-1)/YDIF(J)
  FYM(J)=1.-FY(J)
  ENDDO

  FY(2)=0.
  FYM(2)=1.

  FY(M1)=1.
  FYM(M1)=0.0

    OPEN(UNIT=10,FILE='grid.dat')
    WRITE(10,*) 'VARIABLES="X","Y"'
    WRITE(10,*) 'ZONE I=',L2,',J=',M2,'F=POINT'
     DO J=2,M1
     DO I=2,L1
   	 WRITE(10,*)  XU(I),YV(J)
     ENDDO
    ENDDO
	CLOSE(10)



  END
    
