SUBROUTINE FIRTHCOX(cards, parms, IOARRAY, sort1)
!DEC$ ATTRIBUTES DLLEXPORT :: firthcox
IMPLICIT DOUBLE PRECISION (A-H,O-Z)  

double precision, dimension (15) :: parms
integer N, IP,JCODE,IFIRTH,ISEP,ITER,IMAXIT,IMAXHS
double precision, dimension (int(parms(1))) :: T1, t2
double precision, dimension (int(parms(1)),int(parms(2))) :: X,  bresx
double precision, dimension (int(parms(2)+parms(14))) :: B, B0, FD, absfd, OFFSET,zw1, xx, yy, b0start 
double precision, dimension (int(parms(2)+parms(14)),int(parms(2)+parms(14))) :: SD, VM, WK
integer, dimension (int(parms(1))) :: ibresc, IC
integer, dimension (int(parms(2)+parms(14))) :: IFLAG
double precision, dimension (int(parms(1)),int((2*parms(2)+3+2*(parms(14))))) :: cards
double precision, dimension (int((3+parms(2)+parms(14))),int((parms(2)+parms(14)))) :: IOARRAY
!double precision, dimension (14) :: DER, EREST
logical, dimension (int(parms(2)+parms(14)),int(parms(2)+parms(14))) :: mask
double precision, dimension (int(parms(1)),int(parms(14)+1)) :: ft
integer ngv, ntde
integer, dimension (int(parms(14)+1)) :: ftmap
double precision, dimension (int(parms(1)), int(parms(14)+parms(2))) :: score_weights
integer, dimension (int(parms(1))) :: sort1   ! added by BWJ

INTRINSIC DABS, DSQRT               

!open(unit=6, file="fgcssan.txt")

! ntde = parms(14)
ilike=0
N=int(parms(1))
IP=int(parms(2))
IFIRTH=int(parms(3))
imaxit=int(Parms(4))
imaxhs=int(parms(5))
step=parms(6)
steporig=step
crili=parms(7)
relridge=parms(8)
gconv=parms(9)
if(relridge .lt. 1) relridge=1.

ngv=int(parms(13))
ntde=int(parms(14))
penalty=parms(15)              

parms(10)=-11

offset=ioarray(2,:)
iflag=int(ioarray(1,:))

t1=cards(:,ip+1)
t2=cards(:,ip+2)
ic=int(cards(:,ip+3))
score_weights=cards(:,(ip+4):(2*ip+3+ntde))

x=cards(:,1:ip)
if (ntde .gt. 0) then 
! do j=1,ntde
!  ftmap(j)=ioarray(4,ip+j)
!  write(6,*) ftmap(j)
!  do i=1,n
!   ft(i,j)=cards(i,(2*ip+3+ntde+j))
!   write(6,*) ft(i,j)
!  end do
  ft(:,1:ntde)=cards(:,(2*ip+3+ntde+1):(2*ip+3+ntde*2))
  ftmap(1:ntde)=int(ioarray(4,(ip+1):(ip+ntde)))
!  ftmap(:)=ioarray(4,(ip+1):(ip+ntde))
!  end do
else
 ft=0
 ftmap=0
end if

bresx=x
ibresc=ic


do i=n-1,1,-1
 if (ic(i+1)-1 .gt. -0.0001) then
  if (dabs(t2(i)-t2(i+1)) .lt. 0.0001) then
   score_weights(i,:)=score_weights(i+1,:)
   bresx(i,:)=bresx(i,:)+bresx(i+1,:)
   ibresc(i)=ibresc(i)+ibresc(i+1)
   ibresc(i+1)=0
  end if
 end if
end do

mask=.FALSE.
do j=1,Ip+ntde
 mask(j,j)=.TRUE.
end do

isep=0
XL=0.
xl0=xl-2*crili

b0=0.

where(iflag .eq. 0)
 b0=offset
elsewhere(iflag .eq. 2)
 b0=offset
 iflag=1
endwhere
b0start(:)=b0(:)

isflag=sum(iflag)
!write(6,*) "ISFLAG", isflag
b(:)=b0(:)

ITER=0
iconv=0
JCODE=0


!do while((isflag .gt. 0) .and.  ((iter .gt. 1) .and. (dabs(xl0-xl) .gt. crili)) .and. (iter .lt. imaxit))
do while((iconv .eq. 0) .and. (iter .lt. imaxit))
! if((iter .eq. (imaxit-1)) .and. (imaxit .GT. 10) .and. (step .EQ. steporig)) then
!  step=step/3
!  iter=2
!  b=b0start
! end if
! write(6,*) iter, b
 iter=iter+1
 b0(:)=b(:)
 XL0=XL
 parms(10)=-10
! write(6,*) "Vor 1. LIKE", b
 if (iter .eq. 1) then
  CALL LIKE(N,IP,X,T1,t2,IC,XL,FD,SD,B,JCODE,IFIRTH,ngv,score_weights,bresx,ibresc,ntde,ft,ftmap,penalty,sort1)
  ilike=ilike+1
 end if
! write(6,*) "Nach 1. LIKE"

 parms(10)=-9
 parms(8)= real (JCODE)
 IF (JCODE .GE. 1) RETURN
 IFAIL=0
 wk=-sd
! EPS=.000000000001D0
 CALL INVERT(WK,IP+ntde,IP+ntde,VM,IFAIL)                             
 IF(ITER.EQ.1) then
  parms(12)=xl
  zw=0.
  zw1=matmul(vm, fd)
  zw=dot_product(fd,zw1)
  parms(7)=zw
 end if
 parms(11)=xl
 parms(10)=iter
 parms(9)=isep
 If (ISFLAG.ne.0) then
  IF(IFAIL.ne.0) then
!   "Save" Variance matrix if INVERT failed
 !  WRITE(6,*) 'Inversion failed', ITER,IFIRTH
   ICONV=0
   JCODE=3
   return
  else  
   DO I=1,(IP+ntde)                                                      
    IF (IFLAG(I).EQ.1) then 
     TT=dot_product(vm(I,:),fd(:)*iflag(:))
     IF (DABS(TT).GT.STEP) TT=TT/DABS(TT)*STEP
     B(I)=B(I)+TT
    end if
   end do
!   half step if new log-likelihood is less than old one
   ICONV=0
   IHS=0
   CALL LIKE(N,IP,X,T1,t2,IC,XL,FD,SD,B,JCODE,IFIRTH, ngv, score_weights,bresx,ibresc,ntde,ft,ftmap,penalty,sort1)
   ilike=ilike+1
   wk=-sd
   do while(((XL .le. XL0) .AND. (ITER.ne.1)) .AND. (ihs .le. imaxhs) .and. ((ngv .EQ. IP+ntde) .OR. (ngv .EQ. 0))) 
!   do while((((XL .le. XL0) .or. (sum(abs(fd)) .gt. sumabsfd0)) .AND. (ITER.ne.1)) .AND. (ihs .le. imaxhs) .and. ((ngv .EQ. IP+ntde) .OR. (ngv .EQ. 0))) 
    IHS=IHS+1
    if (abs(relridge-1) .lt. 0.0001) then
     where (iflag .eq. 1)
      b=(b+b0)/2
     end where
    else
     b=b0
     do i=1,(ip+ntde)
      wk(i,i)=wk(i,i)*relridge
     end do
!     EPS=.000000000001D0
     CALL INVERT(WK,IP+ntde,IP+ntde,VM,IFAIL)                             
     DO I=1,(IP+ntde)                                                      
      IF (IFLAG(I).EQ.1) then 
       TT=dot_product(vm(I,:),fd(:)*iflag(:))
       IF (DABS(TT).GT.STEP) TT=TT/DABS(TT)*STEP
        B(I)=B(I)+TT
      end if
     end do
    end if
    CALL LIKE(N,IP,X,T1,t2,IC,XL,FD,SD,B,JCODE,IFIRTH, ngv, score_weights,bresx,ibresc,ntde,ft,ftmap,penalty,sort1)
    ilike=ilike+1
   end do
  end if
 end if
 ICONV=1
 if (isflag .gt. 0) then
  XX=dabs(B-B0)                                                   
  absfd=dabs(fd)*iflag(:) !change 151120GH
  IF(any(XX.GT.CRILI)) ICONV=0
  if(any(absfd .GT. gconv)) ICONV=0
 end if
end do

wk=-sd
!EPS=.000000000001D0


CALL INVERT(WK,IP+ntde,IP+ntde,VM,IFAIL)       


ioarray(3,:)=b(:)
ioarray(4:3+ip+ntde,:)=vm
!do j=1,(ip+ntde)
! ioarray(3,j)=b(j)
! do k=1,(ip+ntde)
!  ioarray(3+j,k)=vm(j,k)
! end do
!end do
!return

!stderr=dsqrt(pack(VM,mask))
!parms(10)=-8



yy=pack(sd,mask)
yy=dabs(yy)
if (any(yy .lt. 0.0001)) isep=1



zw=0.
zw1=matmul(vm, fd*iflag)
zw=dot_product(fd*iflag,zw1)
parms(9)=zw
parms(7)=maxval(dabs(FD)*iflag)
parms(8)=jcode
parms(11)=xl
parms(10)=iter
parms(4)=ilike

!close(unit=6)

RETURN              

end


SUBROUTINE plusone(a)
! this is just a dummy function to check the R-Fortran interface.
double precision a

a=a+1.

return
end



SUBROUTINE INVRT(A,IA)                          


!                                                                       
!...original matrix=a inverse matrix =a (on exit)                                 
!...note that a is changed on exit                                      
!                                                                       
 INTEGER IA,n
 !double precision eps                                             
 double precision, dimension (IA,ia) :: A, B, WK
 INTRINSIC DABS                                                    
                                                                       
 wk=a

 IFAIL=0
 b=a
 N=ia
 
 CALL vert(b, IA, N, WK)
 a=b
    
 RETURN
END  

!SUBROUTINE INVERT(A,IA,N,B,IB,EPS,IFAIL)                          

SUBROUTINE INVERT(A,IA,N,B,IFAIL)                          
!DEC$ ATTRIBUTES DLLEXPORT :: invert


!                                                                       
!...original matrix=a inverse matrix =b                                 
!...note that a is changed on exit                                      
!...eps is a small quantity used to see if matrix singular              
!...ifail on exit ifail=0 ok, ifail=1 matrix nearly singular            
!                                                                       
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 INTEGER IA,N,ifail
 !double precision eps                                             
 double precision, dimension (IA,N) :: A, B, WK
 INTRINSIC DABS                                                    
                                                                       
 wk=a

 IFAIL=0
 b=a

 CALL vert(b, IA, N, WK)

    
 RETURN
END  



SUBROUTINE LIKE(N,IP,X,T1,t2,IC,XL,FD,SD,B,JCODE,IFIRTH, ngv, score_weights,bresx,ibresc,ntde,ft,ftmap,penalty,sort1) 
!DEC$ ATTRIBUTES DLLEXPORT :: like

 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 double precision, dimension (IP+ntde,IP+ntde) :: DINFO, DINFOI, SD, WK
 double precision, dimension (IP+ntde,IP+ntde,IP+ntde) :: dabl, xxxebx
 double precision SEBX, zeitp
 double precision, dimension (IP+ntde) :: XEBX, bresxges
 double precision, dimension (IP+ntde,IP+ntde) :: XXEBX
! double precision, dimension (N+1, IP, Ip, IP) :: XXXEBX
 double precision, dimension (IP+ntde) :: FD, B
 double precision, dimension (N) :: EBX, BX, T1, t2, hh11 !, hh0, hh1, hh2
 double precision hh0, hh1, hh2
 integer, dimension (N) :: IC,ibresc
 double precision, dimension (N,IP) :: X, bresx
 double precision, dimension (N, ip+ntde) :: xges, score_weights
 integer ngv, ntde
 logical, dimension (N) :: maske
 double precision, dimension (N,ntde+1) :: ft
 integer, dimension (ntde+1) :: ftmap
 logical ifastmode
 ! added by BWJ
 integer, dimension (N) :: sort1
 integer indx1, p1
 intrinsic dexp 

 !dlowest=0.000000001
 dlowest=1.0E-44
! write(6,*) "in LIKE"
 XL=0.

 
 ipges=ip+ntde
                                        
 xl=0.
 fd(:)=0.
 sd(:,:)=0.
 dabl(:,:,:)=0.

! Likelihood (XL) is only correct if all or none of the variables is weighted

 xges(:,1:ip)=x

if ((maxval(t1) .lt. minval(t2)) .and. (ntde .eq. 0)) then
  ifastmode=.true.
else
  ifastmode=.false.
end if

if (ifastmode .eqv. .false.) then
 bx=matmul(xges,b)
 ebx=dexp(bx)
 sebx=0.
 xebx=0.
 xxebx=0.
 xxxebx=0.
 indx1=1
 ic_sum=0
 
 do i=N,1,-1
   !look ahead because of Breslow tie correction
   if (i .gt. 1) then   
    if (t2(i) .eq. t2(i-1)) then
     ic_current=0
     ic_sum = ic_sum+ic(i)
    else
     ic_current=ic_sum+ic(i)
     ic_sum=0
    end if
   else
    ic_current=ic_sum+ic(i)
    ic_sum=0
   end if
  sebx = sebx+ebx(i)
  
  do j=1,ipges
   hh0=xges(i,j)*ebx(i)
   xebx(j)=xebx(j)+hh0
   do k=1,ipges
    hh1=hh0*xges(i,k)
	xxebx(j,k)=xxebx(j,k)+hh1
	if (ifirth.eq.1) then
     do l=1,ipges
       hh2=xges(i,l)*hh1
       xxxebx(j,k,l)=xxxebx(j,k,l)+hh2
      end do
     end if
   end do
  end do
	
  ! zeitp=t2(i)-0.00001
  zeitp=t2(i)
  !if (ibresc(i) .ne. 0) then   
  ! bresxges(1:ip)=bresx(i,1:ip)
  if ((ic(i) .ne. 0) .OR. (ic_current .ne. 0)) then
   do while (.true.)   ! subtract the persons with start time > current time
    p1 = sort1(indx1)
    if(t1(p1) .lt. zeitp) then
     exit
	endif
    sebx = sebx-ebx(p1)
	do j=1,ipges
     hh0=xges(p1,j)*ebx(p1)
     xebx(j)=xebx(j)-hh0
     do k=1,ipges
      hh1=hh0*xges(p1,k)
	  xxebx(j,k)=xxebx(j,k)-hh1
	  if (ifirth.eq.1) then
       do l=1,ipges
        hh2=xges(p1,l)*hh1
        xxxebx(j,k,l)=xxxebx(j,k,l)-hh2
       end do
      end if
     end do
    end do
    indx1 = indx1+1;
   end do

   if (sebx .gt. dlowest) then
    dlogsebx=dlog(sebx)
   else
    dlogsebx=dlog(dlowest)
   endif

   if (ngv .eq. ipges) then 
    ! XL=XL+(dot_product(bresxges,b)-ibresc(i)*DLOGSEBX)*score_weights(i,1)
	XL=XL+(bx(i)*ic(i)-ic_current*DLOGSEBX)*score_weights(i,1)
   else
    XL=XL+(bx(i)*ic(i)-ic_current*DLOGSEBX)
	! XL=XL+(dot_product(bresxges,b)-ibresc(i)*DLOGSEBX)
   endif
   
   do j=1,ipges
    ! FD(J)=FD(J)+(bresXges(J)-ibresc(i)*XEBX(J)/SEBX)*score_weights(i,j)
    FD(J)=FD(J)+(xges(i,J)*ic(i)-ic_current*XEBX(J)/SEBX)*score_weights(i,j)	
	do k=1,ipges
     ! SD(J,K)=SD(J,K)-ibresc(i)*((xxebx(j,k)-XEBX(J)/SEBX*XEBX(K))/SEBX)*score_weights(i,j)*score_weights(i,k)
     SD(J,K)=SD(J,K)-ic_current*((xxebx(j,k)-XEBX(J)/SEBX*XEBX(K))/SEBX)*(score_weights(i,j))*(score_weights(i,k))	 
     if (ifirth .ne. 0) then
      ! hh11=xges(:,k)*xges(:,j)*ebx
      DO L=1,IPges
       ! hh2=xges(:,l)*hh11
        DABL(j,k,l)=DABL(j,k,l)-ic_current*((xxxebx(j,k,l)-xxEBX(k,l)*xEBX(j)/SEBx-xEBX(l)*(xxEBX(k,j)  &
         -xEBX(k)*xEBX(j)/SEBx)/SEBx-xEBX(k)*(xxEBX(l,j)-xEBX(l)*xEBX(j)/SEBx)/SEBx)/SEBx)*score_weights(i,j) &
         *score_weights(i,k)*score_weights(i,l)
       ! DABL(j,k,l)=DABL(j,k,l)-ibresc(i)*((xxxebx(j,k,l)-xxEBX(k,l)*xEBX(j)/SEBx-xEBX(l)*(xxEBX(k,j)  &
       ! -xEBX(k)*xEBX(j)/SEBx)/SEBx-xEBX(k)*(xxEBX(l,j)-xEBX(l)*xEBX(j)/SEBx)/SEBx)/SEBx)*score_weights(i,j) &
       ! *score_weights(i,k)*score_weights(i,l)
      end do
     end if
    end do
   end do
   
  end if
 end do
end if

 if (ifastmode .eqv. .true.) then
  xl=0.
  fd=0.
  sd=0.
  xebx=0.
  xxebx=0.
  xxxebx=0.
  sebx=0.
  ic_sum=0
  bx=matmul(xges,b)
  ebx=dexp(bx)
 
  do i=N,1,-1
   if (i .gt. 1) then   !look ahead because of Breslow tie correction
    if (t2(i) .eq. t2(i-1)) then
     ic_current=0
     ic_sum = ic_sum+ic(i)
    else
     ic_current=ic_sum+ic(i)
     ic_sum=0
    end if
   else
    ic_current=ic_sum+ic(i)
    ic_sum=0
   end if
   sebx=sebx+ebx(i)
   do j=1,ipges
    hh0=xges(i,j)*ebx(i)
    xebx(j)=xebx(j)+hh0
    do k=1,ipges
     hh1=hh0*xges(i,k)
     xxebx(j,k)=xxebx(j,k)+hh1
     if (ifirth.eq.1) then
      do l=1,ipges
       hh2=xges(i,l)*hh1
       xxxebx(j,k,l)=xxxebx(j,k,l)+hh2
      end do
     end if
    end do
   end do
   
   if ((ic(i) .ne. 0) .OR. (ic_current .ne. 0)) then
   if (sebx .gt. dlowest) then
    dlogsebx=dlog(sebx)
   else
    dlogsebx=dlog(dlowest)
   endif
   
   if (ngv .eq. ipges) then 
    XL=XL+(bx(i)*ic(i)-ic_current*DLOGSEBX)*score_weights(i,1)
   else
    XL=XL+(bx(i)*ic(i)-ic_current*DLOGSEBX)
   endif
   
   do j=1,ipges
     FD(J)=FD(J)+(xges(i,J)*ic(i)-ic_current*XEBX(J)/SEBX)*score_weights(i,j)                           
     do k=1,ipges
!            SD(J,K)=SD(J,K)-ibresc(i)*((xxebx(j,k)-XEBX(J)/SEBX*XEBX(K))/SEBX)*score_weights(i,j)*score_weights(i,k) 
      SD(J,K)=SD(J,K)-ic_current*((xxebx(j,k)-XEBX(J)/SEBX*XEBX(K))/SEBX)*(score_weights(i,j))*(score_weights(i,k))
      if (ifirth .ne. 0) then
!       hh1=xxebx(j,k)
 !      hh1=xges(:,k)*xges(:,j)*ebx
       DO L=1,IPges
 !       hh2=xges(:,l)*hh1
!        hh2=xges(:,l)*hh1
        DABL(j,k,l)=DABL(j,k,l)-ic_current*((xxxebx(j,k,l)-xxEBX(k,l)*xEBX(j)/SEBx-xEBX(l)*(xxEBX(k,j)  &
         -xEBX(k)*xEBX(j)/SEBx)/SEBx-xEBX(k)*(xxEBX(l,j)-xEBX(l)*xEBX(j)/SEBx)/SEBx)/SEBx)*score_weights(i,j) &
         *score_weights(i,k)*score_weights(i,l)
       end do
      end if
  
     end do
   end do
   end if
  end do
 
 
 
 end if

                    
 
      
    
 if (IFIRTH .NE. 0) then 

  wk(:,:)=-sd(:,:)
  dinfo(:,:)=-sd(:,:)

  IFAIL=0
  CALL INVERT(WK,IPges,IPges,DINFOI,IFAIL)

  DO J=1,IPges
   TRACE=0
   DO  K=1,IPges
    do l=1,IPges
     TRACE=TRACE-DINFOI(k,l)*DABL(j,l,k)
    end do
   end do
   fd(j)=fd(j)+trace*penalty
  end do

 
  DET=DINFO(1,1)
  IF(IPges .NE. 1) then 
   IDETFAIL=0
   DET=0.
!   det=deter(DINFO,IPGES,IPGES)
   det=FindDet(DINFO,IPGES)
!   CALL F03AAF(DINFO,IPges,IPges,DET,WKS,IDETFAIL)
!   IF(IDETFAIL.NE.0) JCODE=2
   IF (DET .LT. dlowest) then 
     DET=dlowest
     JCODE=2
   end if
  end if
  XL=XL+penalty*dlog(DET)
 end if

 RETURN
END                                                               





SUBROUTINE PLCOMP(CARDS, PARMS, IOARRAY, sort1)
!DEC$ ATTRIBUTES DLLEXPORT :: plcomp

IMPLICIT DOUBLE PRECISION (A-H,O-Z)
INTRINSIC DABS, DSQRT, DSIGN                                                    
double precision, dimension (15) :: parms
double precision, dimension (15) :: parmsfc
double precision, dimension (int(parms(1))) :: T1,t2
double precision, dimension (int(parms(1)),int(parms(2))) :: X, bresx
integer, dimension (int(parms(1))) :: IC, ibresc
double precision, dimension (int(parms(2)+parms(14))) :: B, B0, FD, OFFSET, PVALUE, XE,DELTA,XGRAD, BSAVE
double precision, dimension (int(parms(2)+parms(14)),int(parms(2)+parms(14))) :: SD, VM, WK, XHESS, XVM
double precision, dimension (int(parms(2)+parms(14)),2) :: CI
integer, dimension (int(parms(2)+parms(14))) :: IFLAG
double precision, dimension (int(parms(1)),int(2*parms(2)+2*parms(14)+3)) :: cards
double precision, dimension (9,int(parms(2)+parms(14))) :: IOARRAY
double precision, dimension (int(3+parms(2)+parms(14)),int(parms(2)+parms(14))) :: IOAFC
double precision, dimension (int(parms(1)),int(parms(14)+1)) :: ft
integer, dimension (int(parms(14)+1)) :: ftmap
double precision, dimension (int(parms(1)), int(parms(14)+parms(2))) :: score_weights
integer, dimension (int(parms(1))) :: sort1   ! added by BWJ

!open(unit=6, file="fgcsspl.txt")

N=int(parms(1))
IP=int(parms(2))
IFIRTH=int(parms(3))
imaxit=int(Parms(4))
imaxhs=int(parms(5))
step=parms(6)
crili=parms(7)
CHI=parms(8)
gconv=PARMS(9)
PARMS(9)=0
ngv=int(parms(13))
ntde=int(parms(14))
penalty=parms(15)

offset(:)=ioarray(2,:)
iflag(:)=int(ioarray(1,:))
b(:)=ioarray(3,:)
b0(:)=ioarray(3,:)

!   do j=1,ip,1
!    offset(j)=IOARRAY(2,j)
!    iflag(j)=IOARRAY(1,j)
!    b(j)=IOARRAY(3,j)
!    b0(j)=IOARRAY(3,j)
!   end do

t1=cards(:,ip+1)
t2=cards(:,ip+2)
ic=int(cards(:,ip+3))


x=cards(:,1:ip)
if (ntde .gt. 0) then 
  ft(:,1:ntde)=cards(:,(2*ip+3+ntde+1):(2*ip+3+ntde*2))
  ftmap(1:ntde)=int(ioarray(4,(ip+1):(ip+ntde)))
else
 ft=0
 ftmap=0
end if


bresx=x
ibresc=ic

do i=n-1,1,-1
 if (ic(i+1)-1 .gt. -0.0001) then
  if (dabs(t2(i)-t2(i+1)) .lt. 0.0001) then
   score_weights(i,:)=score_weights(i+1,:)
   bresx(i,:)=bresx(i,:)+bresx(i+1,:)
   ibresc(i)=ibresc(i)+ibresc(i+1)
   ibresc(i+1)=0
  end if
 end if
end do


score_weights=cards(:,(ip+4):(2*ip+3+ntde))



!   do i=1,n,1
!    t(i)=cards(i,ip+1)
!    ic(i)=cards(i,ip+2)
!    do j=1,ip,1
!     x(i,j)=cards(i,j)
!    end do
!   end do



!   XL ... log likelihood
!   FD ... First Derivative of log likelihood (=score vector)
!   SD ... Second Derivative of log likelihood
!   VM ... Variance Matrix ( = -SD**(-1) )
!   DINFO .Fisher Information (= -SD )
!   DINFOI.Inverse Fisher Information ( = VM )
!   SDI ...Inverse seconde derivative 

!   BX  ...Xb (n x 1)
!   EBX ...exp(Xb) (n x 1)
!   XEBX ..X'exp(Xb) (IP x 1)
!   XXEBX .sum(i)sum(j)sum(k)(X(i,k)*X(i,j)*exp(X(i,)*b))

!   XE  ...Unit vector
!     CI  ...confidence interval

    
!   Assuming that B maximizes the penalized likelihood

!write(6,*) "vor 1. LIKE"
CALL LIKE(N,IP,X,T1,t2,IC,XL,FD,SD,B,JCODE,IFIRTH,ngv,score_weights,bresx,ibresc, ntde,ft,ftmap,penalty,sort1) 
!write(6,*) "nach 1. Like"
wk(:,:)=-sd(:,:)
!EPS=.000000000001D0
CALL INVERT(WK,IP+ntde,IP+ntde,VM,IFAIL)                             
XLMAX=XL
IFAIL=0
!CONFLEV=1.-ALPHA
DCRILI=CRILI
DF=1.
!CHI=G01FCF(CONFLEV,DF,IFAIL)
XL0=XLMAX-0.5*CHI

bsave(:)=b(:)
xgrad(:)=fd(:)
xvm(:,:)=vm(:,:)
xhess(:,:)=sd(:,:)


!   DO 100 K=1,IP
do k=1,ip+ntde
!DO 101 L=1,2
 IF (iflag(k) .EQ. 1) then
     do L=1,2

    !   L=1 ... lower limit
    !   L=2 ... upper limit

      xe(:)=0

      XE(K)=1

      XSIGN=L*2-3
      XLAMBDA=0
    !   Initialization of Likelihood, first and second der., beta
      vm(:,:)=xvm(:,:)
      sd(:,:)=xhess(:,:)
      fd(:)=xgrad(:) * iflag(:) ! change 151120GH
      b(:)=bsave(:)

        
      ICONV=0
      ITER=0

    !   DO 140 WHILE (ICONV .EQ. 0)
      DO WHILE (ICONV .EQ. 0)
       ITER=ITER+1  
       DO K2=1,IP+ntde
        DELTA(K2)=0

        DO K3=1,IP+ntde
         DELTA(K2)=DELTA(K2)+(FD(K3)+XLAMBDA*XE(K3))*VM(K2,K3)
        end do
        IF (DABS(DELTA(K2)) .GT. step) THEN
         DELTA(K2)=DSIGN(step,DELTA(K2))
        ENDIF
        B(K2)=B(K2)+DELTA(K2)
       end do

       CALL LIKE(N,IP,X,T1,t2,IC,XL,FD,SD,B,JCODE,IFIRTH,ngv,score_weights,bresx,ibresc,ntde,ft,ftmap,penalty,sort1) 
       wk(:,:)=-sd(:,:)
       fd(:) = fd(:) * iflag(:) !change 151120GH
!       EPS=.000000000001D0
       CALL INVERT(WK,IP+ntde,IP+ntde,VM,IFAIL)                             


       GVG=0
       DO K2=1,IP+ntde
        DO K3=1,IP+ntde
!         GVG=GVG+FD(K2)*FD(K3)*VM(K2,K3)   !commented 151120GH
         GVG=GVG+FD(K2)*FD(K3)*VM(K2,K3)*iflag(K2)*iflag(K3)   !151120GH
        end do
       end do

       XLAMBDA=XSIGN*(DSQRT(-(2*(XL0-XL-0.5*GVG)/VM(K,K))))


       AXLDIFF=DABS(XL0-XL) 
       IF (AXLDIFF <= DCRILI) THEN 
        ICONV=1
       ENDIF

       C2=0

       DO K2=1,IP+ntde
        DO K3=1,IP+ntde
!         C2=C2+(FD(K2)+XLAMBDA*XE(K2))*(VM(K2,K3))*(FD(K3)+XLAMBDA*XE(K3))           !commented 151120GH
         C2=C2+(FD(K2)+XLAMBDA*XE(K2))*(VM(K2,K3))*(FD(K3)+XLAMBDA*XE(K3))*iflag(K3)*iflag(K2) !151120GH
        end do
       end do
        
       IF (C2 .GT. DCRILI) THEN 
        ICONV=0
       ENDIF

       IF (ITER .EQ. IMAXIT) THEN 
        ICONV=1
       ENDIF

      end do

      CI(K,L)=B(K)

      IOARRAY(6+L,K)=ITER
        
     end do
 endif
end do      

!   penalized LR p-values for H0:b(k)=0
    
!write (6,*) "vor p-value berechnung"

DO K2=1,IP+ntde
 if (IFLAG(K2) .EQ. 1) then
     OSSAVE=OFFSET(K2)
     IFLSAVE=IFLAG(K2)
     OFFSET(K2)=0
!     IFLAG(K2)=0
     ITER=0
    !   CALL FIRTHCOX(N,IP,X,T,IC,B,B0,XL,FD,SD,VM,BX,           
    !     +JCODE,IFIRTH,CRILI,ISEP,ITER,IMAXIT,IMAXHS,STEP,IFLAG,
    !     +OFFSET)

     parmsfc(1)=n
     parmsfc(2)=ip
     parmsfc(3)=ifirth
     parmsfc(4)=imaxit
     parmsfc(5)=imaxhs
     parmsfc(6)=step
     parmsfc(7)=crili
     parmsfc(8)=1
     parmsfc(9)=gconv
     parmsfc(10)=0
     parmsfc(11)=0
     parmsfc(12)=0
     parmsfc(13)=ngv
     parmsfc(14)=ntde
     parmsfc(15)=penalty

     IOAFC(2,:)=bsave(:)*IFLAG(:)   !151120GH, starting values for p-value computation
        
     do jjj=1,IP+ntde,1
    !  IOAFC(1,jjj)=1               !commented 151120GH
       IOAFC(1,jjj)=IFLAG(jjj)      !changed 151120GH
!       IOAFC(2,jjj)=0.
     end do
    !   IFLAG for estimation set to 0
     IOAFC(1,K2)=0.
    !   offset value 0
     IOAFC(2,K2)=0.
     if (ntde .gt. 0) then
      do jjj=1, ntde
       ioafc(4,ip+jjj)=ftmap(jjj)
      end do
     end if
     CALL FIRTHCOX(CARDS, PARMSFC, IOAFC, sort1)
     IF (PARMSFC(8).GE.1) PARMS(9)=2
    ! write(6,*) "Var: " , k2
    ! write(6,*) "Code: ", parmsfc
    ! write(6,*) "likelihood: ", xl, xlmax
     xl=parmsfc(11)
     CHI=2*(XLMAX-XL)       
     IFAIL=0
     DF=1
    ! PV=G01ECF('U',CHI,DF,IFAIL)
     PVALUE(K2)=CHI                     ! it is not a p-value, just the chi-square
     OFFSET(K2)=OSSAVE
!     IFLAG(K2)=IFLSAVE
    ioarray(9,k2)=parmsfc(10)
 else
    pvalue(k2)=0.
    ioarray(9,k2)=0
 endif
end do



vm(:,:)=xvm(:,:)
sd(:,:)=xhess(:,:)
fd(:)=xgrad(:)
b(:)=bsave(:)
    
do j=1,ip+ntde
 do l=1,2
  ioarray(3+l,j)=ci(j,l)
 end do
 ioarray(6,j)=pvalue(j)
end do
    
!   JCODES (PARMS(9)): 2=problem in FIRTHCOX-CALL


!close(unit=6)

RETURN
END   


!

! Code converted using TO_F90 by Alan Miller
! Date: 2006-04-13  Time: 10:10:03

!      ________________________________________________________
!     |                                                        |
!     |     FACTOR A GENERAL MATRIX WITH PARTIAL PIVOTING      |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         A     --ARRAY CONTAINING MATRIX                |
!     |                 (LENGTH AT LEAST 3 + N(N+1))           |
!     |                                                        |
!     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
!     |                                                        |
!     |         N     --DIMENSION OF MATRIX STORED IN A        |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         A     --FACTORED MATRIX                        |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS                              |
!     |    PACKAGE SUBROUTINES: PACK                           |
!     |________________________________________________________|

SUBROUTINE fact(ain,a,n)

!INTEGER, INTENT(IN)                      :: la
INTEGER, INTENT(IN)                      :: n
double precision, INTENT(in OUT)                     :: a(3+n*(n+1))
double precision, intent(in)              ::  ain(n,n)
double precision :: r,s,t
INTEGER :: e,f,g,h,i,j,k,l, m, o,p

intrinsic dabs

!IF ( la > n ) CALL pack(a,la,n)
!nicht notwendig da la immer = N

do i=1,n
 do j=1,n
  a((i-1)*n+j)=ain(i,j)
 end do
end do


r = 0.
o = n + 1
p = o + 1
l = 5 + n*p
i = -n - 3
!     ---------------------------------------------
!     |*** INSERT PIVOT ROW AND COMPUTE 1-NORM ***|
!     ---------------------------------------------
10    l = l - o
IF ( l == 4 ) GO TO 30
s = 0.
DO  k = 1,n
  j = l - k
  t = a(i+j)
  a(j) = t
  s = s + DABS(t)
END DO
IF ( r < s ) r = s
i = i + 1
GO TO 10
30    a(1) = 1230
a(2) = n
a(3) = r
i = 5 - p
k = 1
40    i = i + p
IF ( k == n ) GO TO 110
e = n - k
m = i + 1
h = i
l = i + e
!     ---------------------------------------
!     |*** FIND PIVOT AND START ROW SWAP ***|
!     ---------------------------------------
DO  j = m,l
  IF ( DABS(a(j)) > DABS(a(h)) ) h = j
END DO
g = h - i
j = i - k
a(j) = g + k
t = a(h)
a(h) = a(i)
a(i) = t
k = k + 1
IF ( t == 0. ) GO TO 100
!     -----------------------------
!     |*** COMPUTE MULTIPLIERS ***|
!     -----------------------------
DO  j = m,l
  a(j) = a(j)/t
END DO
f = i + e*o
70    j = k + l
h = j + g
t = a(h)
a(h) = a(j)
a(j) = t
l = e + j
IF ( t == 0. ) GO TO 90
h = i - j
!     ------------------------------
!     |*** ELIMINATE BY COLUMNS ***|
!     ------------------------------
m = j + 1
DO  j = m,l
  a(j) = a(j) - t*a(j+h)
END DO
90    IF ( l < f ) GO TO 70
GO TO 40
100   a(1) = -1230
GO TO 40
110   IF ( a(i) == 0. ) a(1) = -1230
RETURN
END SUBROUTINE fact





!

! Code converted using TO_F90 by Alan Miller
! Date: 2006-04-13  Time: 10:09:46

!      ________________________________________________________
!     |                                                        |
!     |                INVERT A GENERAL MATRIX                 |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         V     --ARRAY CONTAINING MATRIX                |
!     |                                                        |
!     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
!     |                                                        |
!     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
!     |                                                        |
!     |         W     --INTEGER WORK ARRAY WITH AT LEAST N-1   |
!     |                      ELEMENTS                          |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         V     --INVERSE                                |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS                              |
!     |________________________________________________________|

SUBROUTINE vert(v,lv,n,w)

INTEGER, INTENT(IN OUT)                  :: lv
INTEGER, INTENT(IN)                      :: n
double precision, INTENT(IN OUT)                     :: v(lv,N)
double precision, INTENT(OUT)                     :: w(N)
double precision :: s,t
INTEGER :: i,j,k,l,m, p

!Anm GH bei den Dimensionen die Indizes verändert, waren v(lv,1) und w(1) vorher

intrinsic dabs
!GH 161221
k = 1

IF ( n == 1 ) GO TO 110
l = 0
m = 1
10    IF ( l == n ) GO TO 90
k = l
l = m
m = m + 1
!     ---------------------------------------
!     |*** FIND PIVOT AND START ROW SWAP ***|
!     ---------------------------------------
p = l
IF ( m > n ) GO TO 30
s = DABS(v(l,l))
DO  i = m,n
  t = DABS(v(i,l))
  IF ( t <= s ) CYCLE
  p = i
  s = t
END DO
w(l) = p
30    s = v(p,l)
v(p,l) = v(l,l)
IF ( s == 0. ) GO TO 120
!     -----------------------------
!     |*** COMPUTE MULTIPLIERS ***|
!     -----------------------------
v(l,l) = -1.
s = 1./s
DO  i = 1,n
  v(i,l) = -s*v(i,l)
END DO
j = l
50    j = j + 1
IF ( j > n ) j = 1
IF ( j == l ) GO TO 10
t = v(p,j)
v(p,j) = v(l,j)
v(l,j) = t
IF ( t == 0. ) GO TO 50
!     ------------------------------
!     |*** ELIMINATE BY COLUMNS ***|
!     ------------------------------
IF ( k == 0 ) GO TO 70
DO  i = 1,k
  v(i,j) = v(i,j) + t*v(i,l)
END DO
70    v(l,j) = s*t
IF ( m > n ) GO TO 50
DO  i = m,n
  v(i,j) = v(i,j) + t*v(i,l)
END DO
GO TO 50
!     -----------------------
!     |*** PIVOT COLUMNS ***|
!     -----------------------
90    l = int(w(k))
DO  i = 1,n
  t = v(i,l)
  v(i,l) = v(i,k)
  v(i,k) = t
END DO
k = k - 1
IF ( k > 0 ) GO TO 90
RETURN
110   IF ( v(1,1) == 0. ) GO TO 120
v(1,1) = 1./v(1,1)
RETURN
120   continue
!WRITE(6,*) 'ERROR: MATRIX HAS NO INVERSE'
!STOP
END SUBROUTINE vert




!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
DOUBLE PRECISION FUNCTION FindDet(matrix, n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, DIMENSION(n,n) :: matrix
    DOUBLE PRECISION :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
    
END FUNCTION FindDet
