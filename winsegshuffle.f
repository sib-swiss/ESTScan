*       Program winsegshuffle
*       $Id: winsegshuffle.f,v 1.2 2005/03/10 15:08:04 c4chris Exp $
*----------------------------------------------------------------------*     
*       Function: generates a window-segment shuffled database
*		  from a source sequence database 
*       Author:   Philipp Bucher
*----------------------------------------------------------------------*     

	Parameter        (NERR=   0)

        Character*01      CSEQ(11000000)
        Character*01      CSEG(100)
	Integer           IWIN(100000)

	Character*132     RCIN


	Character*64      CHHA

* Read command line 

        NTOT=10000000
        NLEN=200
	NSEG=10
	NWIN=20
        LINL=70
	IDUM=-100

        Call Repar
     *    (NTOT,NLEN,NSEG,NWIN,IDUM,IRC)
        If(IRC.NE.0) then
           Write(NERR,'(
     *      ''Usage: winsegshuffle '',
     *      ''[db_size [seq_length [seg_length '',
     *      ''[window [iran-seed]]]]]''
     *        )')
           Stop
        End if

C       Print *,NTOT
C       Print *,NLEN
C       Print *,NSEG
C       Print *,NWIN
C       Print *,IDUM

* Read input sequences

        K1=0
    1   Read(5,'(A)',End=40) RCIN
	If(RCIN(1:1).EQ.'>') then
	   If(K1.GE.NTOT) Go to  50
        Else 
	   L=Lblnk(RCIN)
           Read(RCIN,'(132A)')(CSEQ(ii1),ii1=K1+1,K1+L)
           K1=K1+L
        End if  	
	Go to   1

   40   NTOT = K1 - L

* Do segment window shuffling

   50   Continue

        Do I1=1,NTOT,NSEG*NWIN

           Call Permut(IWIN,NWIN,IDUM)
*          Write(6,'(20I3)')(IWIN(ii1),ii1=1,20)
	         K3=0
	   Do I2=1,NWIN
	         N3=I1+(IWIN(I2)-1)*NSEG
	      Do I3=1,NSEG
                 K3=K3+1
	         N3=N3+1
		 CSEG(K3)=CSEQ(N3)
              End do 
           End do
	      K2=0
	   Do I2=I1+1,I1+NSEG*NWIN
	      K2=K2+1
	      CSEQ(I2)=CSEG(K2)
           End do
        End do

* Print sequence:

           K1=0
        Do I1=1,NTOT,NLEN
	   K1=K1+1
           Write(CHHA,'(''>'',2I5)') NLEN,K1
	   Do I2=2,11
	     If(CHHA(I2:I2).EQ.' ') CHHA(I2:I2)='0'
             End do 
           CHHA=CHHA(1:11) // ' ..'

           Write(6,'(64A)')(CHHA(ii1:ii1),ii1=1,Lblnk(CHHA))
           Write(6,'((70A))')(CSEQ(ii1),ii1=I1+1,I1+NLEN)

	End do

  100   Stop
  900   Go to 100 

        End
*----------------------------------------------------------------------*
        FUNCTION RAN2(IDUM)
        PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112E-6)
        DIMENSION IR(97)
        DATA IFF /0/
        IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
          IFF=1
          IDUM=MOD(IC-IDUM,M)
          DO 11 J=1,97
            IDUM=MOD(IA*IDUM+IC,M)
            IR(J)=IDUM
   11     CONTINUE
          IDUM=MOD(IA*IDUM+IC,M)
          IY=IDUM
        ENDIF
        J=1+(97*IY)/M
        IF(J.GT.97.OR.J.LT.1) PAUSE
        IY=IR(J)
        RAN2=IY*RM
        IDUM=MOD(IA*IDUM+IC,M)
        IR(J)=IDUM
        RETURN
        END
*----------------------------------------------------------------------*     
        Subroutine Permut(IR,NR,IDUM)

        Integer IR(*)
           Do  10 I1=1,NR 
              IR(I1)=I1
   10      Continue

           Do  20 I1=NR,2,-1
              J1=IR(I1)
              RS=RAN2(IDUM)
              K1=INT(RS*I1)+1
              IR(I1)=IR(K1)
              IR(K1)=J1
   20      Continue

        Return
        End
*----------------------------------------------------------------------*     
        Function          Lblnk(string) 
        Character*(*)     string 

        L=Len(string)

        Do   9 I1=L,1,-1
           If(STRING(I1:I1).NE.' ') go to  10
    9   Continue
   10   Lblnk=I1

        Return
        End
*----------------------------------------------------------------------*     
        Subroutine Repar
     *    (NTOT,NLEN,NSEG,NWIN,IDUM,IRC)

	Character*64      CARG

	IRC=0

        N1=Iargc()
	Do I1=1,N1
           Call GetArg(I1,CARG)
	   If(CARG(1:1).NE.'-') then
              If(I1.EQ.1) Read(CARG,*,Err=900) NTOT
              If(I1.EQ.2) Read(CARG,*,Err=900) NLEN
              If(I1.EQ.3) Read(CARG,*,Err=900) NSEG
              If(I1.EQ.4) Read(CARG,*,Err=900) NWIN
              If(I1.EQ.5) Read(CARG,*,Err=900) IDUM
           End if
        End do 

        If(IDUM.GT.0) IDUM=-1-IDUM

  100   Return 
  900   IRC=1
        Go to 100
        End
*----------------------------------------------------------------------*
