    subroutine butter(frq,fs,iband,l,bn,an)
!----------------------------------------------------------------------
! Routine butter:To design LP,HP,BP,BS digital Butterworth filter.
! All array has been specified internally,so user only need to
! input f1,f2,f3,f4,fs(in hz), IBAND (to specify the type of to design)
! and l (the order of filter). This program output Hk(z)=Bk(z)/Ak(z),k=1,2,...,
! ksection.
! Input parameters:
!                               l: the order of filter
!                           f1,f2: the edge frequency desired(in Hz):
!      For LP IIR Filter(iband=1): f1=frq(1),f2=0
!      For HP IIR Filter(iband=2): f1=frq(1),f2=0
!      For BP IIR Filter(iband=3): f1=frq(1),f2=frq(2)
!      For BS IIR Filter(iband=4): f1=frq(1),f2=frq(2)
!                       fs: sampling frequnecy in Hz.
!  Output parameters: 
!                              bn: the coefficient of numerator
!                              an: the coefficient of denominator
!   Note: the screen output will demostrate your IIR system function
!         in cascade form.
!                               modified from Hu Guangshu in chapter 7
!----------------------------------------------------------------------
    implicit none
    real*8,allocatable :: frq(:)
    real*8 :: fs
    integer :: l,iband
    
    real*8,allocatable :: bn(:),an(:)
    real*8,allocatable :: pb(:),pa(:)
    real*8,allocatable :: tmpb(:,:),tmpa(:,:)
    
    real*8 :: b(0:4),a(0:4),d(0:4),c(0:4)
    integer :: ln=2
    
    integer :: lo,lh,ksection,ierror
    integer :: i,k,k1,l1,l2
    real*8 :: f1,f2
    real*8 :: fln,flh,t1
    
    !print *, f1,f2,f3,f4,fs
    
    if(iband == 1) then
!------------------  for low pass ------------------------
        if (size(frq) /= 1) then
            print *, 'the number of cutoff frequency for LP is wrong'
            return
        end if
        f1 = frq(1)
        f2 = 0.0
        !f1 = 0.0
    end if
    
    if (iband == 2) then
!------------------  for high pass ------------------------  
        if (size(frq) /= 1) then
            print *, 'the number of cutoff frequency for HP is wrong'
            return
        end if
        f1 = frq(1)
        f2 = 0.0    
    end if
    if (iband == 3) then
!------------------  for band pass ------------------------
        if (size(frq) /= 2) then
            print *, 'the number of cutoff frequency for BP is wrong'
            return
        end if
        f1 = frq(1)
        f2 = frq(2)
    end if
    if (iband == 4) then
!------------------  for band stop ------------------------
        if (size(frq) /= 2) then
            print *, 'the number of cutoff frequency for BS is wrong'
            return
        end if
        f1 = frq(1)
        f2 = frq(2)
    end if
    
    !print *, f1,f2,fs
    
    if((iband == 1) .or. (iband == 2)) lo = l
    if((iband == 3) .or. (iband == 4)) lo = l*2
    allocate(bn(0:lo))
    allocate(an(0:lo))
    allocate(pb(0:lo))
    allocate(pa(0:lo))
    
    fln = f1/fs
    flh = f2/fs
    if (l == 1) then
        ksection = 1
    end if
    t1 = l/2.0
    lh = int(l/2)
    if (t1 /= float(lh)) ksection = (l+1)/2
    if (t1 == float(lh)) ksection = l/2
    if ((iband == 3) .or. (iband == 4)) ln = 4
    
    allocate(tmpb(1:ksection,0:ln))
    allocate(tmpa(1:ksection,0:ln))
    
    do k1=1,ksection
        call butwcf(l,k1,ln,d,c,ierror)
        if (ierror /= 0) then
            print *, 'stop at routine BUTWCF, ierror = ',ierror
            return
        end if
!        write(*,*)'--------------- for ksection=',k1,'------------'
!        write(*,*)'    Analog low-pass filter hk(s)=dk(s)/ck(s)'
!        do k=0,ln
!            write(*,70)k,d(k),k,c(k)
!        end do
!70      format('     d(',i2,')=',f10.4,'        c(',i2,')=',f10.4)
        call aftodf(d,c,ln,iband,fln,flh,b,a,ierror)
        if (ierror /= 0) then
            print *,'stop at routine AFTODF, ierror = ',ierror
            return
        end if
        a(0) = 1.0
!        write(*,*)'    Digital low-pass filter hk(z)=dk(z)/ck(z)'
        do k=0,ln
            tmpb(k1,k) = b(k)
            tmpa(k1,k) = a(k)
!            write(*,90)k,b(k),k,a(k)
        end do
!90      format('     b(',i2,')=',f10.4,'        a(',i2,')=',f10.4)
        
    end do
    
!---------------collect the coefficients of each section----------
    if (t1 /= float(lh)) then
        if((iband == 3) .or. (iband==4)) l1 = 2
        if((iband == 1) .or. (iband==2)) l1 = 1
    else
        l1 = ln
    end if
    do k=0,l1
        bn(k) = tmpb(ksection,k)
        an(k) = tmpa(ksection,k)
        pb(k) = tmpb(ksection,k)
        pa(k) = tmpa(ksection,k)
    end do
    l2 = ln
    do k1=1,ksection-1
        do k=0,ln
            b(k) = tmpb(k1,k)
            a(k) = tmpa(k1,k)
        end do
        call collect(pb,b,l1,l2,bn)
        call collect(pa,a,l1,l2,an)
        l1 = l1+ln
        l2 = ln
        do k=0,l1
            pb(k) = bn(k)
            pa(k) = an(k)
        end do
    end do
    
    deallocate(tmpb)        
    deallocate(tmpa)
    deallocate(pb)
    deallocate(pa)
    
    l = lo
    
    return
    end subroutine butter
    
    subroutine collect(p1,p2,l1,l2,pn)
	implicit none
	integer :: l1,l2
	real*8 :: p1(0:l1),p2(0:l2),pn(0:l1+l2)
 
	integer :: i,j
	
	do i=0,l1+l2
		pn(i)=0.0
	end do
 
	do i=0,l2
		do j=0,l1
			pn(j+i)=pn(j+i)+p2(i)*p1(j)
		end do
	end do
 
	return
    end subroutine collect
    
    subroutine biline(work,d,c,ln,b,a,ierror)
!----------------------------------------------------------------------
! Routine BILINE: To convert analog H(S) to digital H(Z) via bilinear
!                 transformation. H(S)=D(S)/C(S), H(Z)=B(Z)/A(Z)
! LN specifies the length of the coefficient arrays and filter order L
! is computed internally.   WORK is an internal array (2D)
!   IF  IERROR=0:    no errors detected in transformation
!             =1:    all zero transfer function
!             =2:    invalid transfer function; y(k) coef=0
!       From Ref. [5] of Chapter 2 .       in chapter 7
!----------------------------------------------------------------------
    implicit none
    real*8 :: work(0:4,0:4),d(0:4),c(0:4),b(0:4),a(0:4)
    integer :: ln,ierror
    
    integer :: iflag
    integer :: i,j,l
    real*8 :: tmp,atmp,scale
    
    iflag = 1
    do i=ln,0,-1
        if ((c(i) /= 0.0) .or. (d(i) /= 0.0)) then
            iflag = 0
            exit
        end if
    end do
    if (iflag == 1) then
        ierror = 1
        return
    end if
    l = i
    
    do j=0,l
        work(0,j) = 1.0
    end do
    tmp = 1.0
    do i=1,l
        tmp = tmp*float(l-i+1)/float(i)
        work(i,0) = tmp
    end do
    do i=1,l
        do j=1,l
            work(i,j) = work(i,j-1)-work(i-1,j)-work(i-1,j-1)
        end do
    end do
    do i=l,0,-1
        b(i) = 0.0
        atmp = 0.0
        do j=0,l
            b(i) = b(i)+work(i,j)*d(j)
            atmp = atmp+work(i,j)*c(j)
        end do
        scale = atmp
        if(i /= 0) a(i) = atmp
    end do
    ierror = 2
    if(scale == 0.0) return
    b(0) = b(0)/scale
    do i=1,l
        b(i) = b(i)/scale
        a(i) = a(i)/scale
    end do
    do i=l+1,ln
        b(i) = 0.0
        a(i) = 0.0
    end do
    ierror = 0
    return
    end subroutine biline
    
    subroutine aftodf(d,c,ln,iband,fln,fhn,b,a,ierror)
!----------------------------------------------------------------------
! Routine AFTODF: To convert normalized LP analog H(s) to digital H(z).
!   H(s)=D(s)/C(s),H(z)=B(z)/A(z).Filter order l is computed internally.
!   LN specifies coefficient array size. WORK(0:LN,0:LN) is a work array.
!   IF   IBAND=1:    lowpass   fln=normalized cutoff frequency
!             =2:    highpass  fln=normalized cutoff frequency
!             =3:    bandpass  fln=low  cutoff frequency
!                              fhn=high cutoff frequency
!             =4:    bandstop  fln=low  cutoff frequency
!                              fhn=high cutoff frequency
!   IF  IERROR=0:    no errors detected
!              1:    all zero transfer function
!              2:    biline: invalid transfer function
!              3:    filter order exceeds array size
!              4:    invalid filter type parameter (IBAND)
!              5:    invalid cutoff frequency
!       From Ref. [5] of Chapter 2 .      in chapter 7
!-----------------------------------------------------------------------
    implicit none
    real*8 :: work(0:4,0:4),d(0:4),c(0:4),b(0:4),a(0:4)
    integer :: ln,iband,ierror
    real*8 :: fln,fhn
    
    real*8 :: pi,w,w1,w2,w02
    real*8 :: tmpd,tmpc,tmp
    integer :: i,k,ll,mm,ls
    integer :: iflag,m,l
    
    real*8 :: spbfct
    
    pi = 4.0*atan(1.0)
    ierror = 0
    if((iband < 1) .or. (iband >4)) ierror = 4
    if((fln <= 0.0) .or. (fln > 0.5)) ierror = 5
    if((iband >= 3) .and. (fln >= fhn)) ierror = 5
    if((iband >= 3) .and. (fhn > 0.5)) ierror = 5
    if(ierror /= 0) return
    
    iflag = 1
    do i=ln,0,-1
        if ((c(i) /= 0.0) .or. (d(i) /= 0.0)) then
            iflag = 0
            exit
        end if
    end do
    if (iflag == 1) then
        ierror = 1
        return
    end if
    
    m = i    
    w1 = tan(pi*fln)
    l = m
    if (iband > 2) then
        l = 2*m
        w2 = tan(pi*fhn)
        w = w2-w1
        w02 = w1*w2
    end if
    ierror = 3
    if(l > ln) return
    
    if ((iband == 2) .or. (iband == 4)) then
!-------- substitution of 1/s to generate highpass (hp,bs) ------------
        do mm=0,m/2
            tmp = d(mm)
            d(mm) = d(m-mm)
            d(m-mm) = tmp
            tmp = c(mm)
            c(mm) = c(m-mm)
            c(m-mm) = tmp
        end do
    end if
    
    if ((iband == 1) .or. (iband == 2)) then
!-------- scaling s/w1 for lowpass, highpass ---------------------------
        do mm=0,m
            d(mm) = d(mm)/(w1**mm)
            c(mm) = c(mm)/(w1**mm)
        end do
    end if
            
    if (iband == 3 .or. iband == 4) then
!-------- substitution of (s**2+w0**2)/(w*s)  bandpass,bandstop -------
        do ll=0,l
            work(ll,1) = 0.0
            work(ll,2) = 0.0
        end do
        do mm=0,m
            tmpd = d(mm)*(w**(m-mm))
            tmpc = c(mm)*(w**(m-mm))
            do k=0,mm
                ls = m+mm-2*k-1
                !print *, spbfct(mm-k,mm-k),mm-k
                tmp = spbfct(mm,mm)/(spbfct(k,k)*spbfct(mm-k,mm-k))
                work(ls+1,1) = work(ls+1,1)+tmpd*(w02**k)*tmp
                work(ls+1,2) = work(ls+1,2)+tmpc*(w02**k)*tmp
            end do
        end do
        do ll=0,l
            d(ll) = work(ll,1)
            c(ll) = work(ll,2)
        end do
    end if
    
!---------- substitute (z-1)/(z+1) ------------------------------------
    call biline(work,d,c,ln,b,a,ierror)
    if (ierror /= 0) then
        print *, 'stop at routine biline, ierror = ',ierror
    end if
    
    return
    end subroutine aftodf
    
    function spbfct(i1,i2)
!-------- generates (i1)!/(i1-i2)!=i1*(i1-1)*...*(i1-i2+1). -----------
!-------- note: 0!=1 and spbfct(i,i)=spbfct(i,i-1)=i!.      -----------
    integer :: i1,i2
    real*8 :: spbfct
    
    spbfct = 0.0
    if((i1 < 0) .or. (i2 < 0) .or. (i2 > i1)) return
    spbfct = 1.0
    if(i2 == 0) return
    do i=i1,i1-i2+1,-1
        spbfct = spbfct*i
    end do
    
    return
    end function spbfct
    
    subroutine butwcf(l,k,ln,d,c,ierror)
!----------------------------------------------------------------------
!    routine BUTWCF: To design low-pass Butterworth analog filter:
!                    H(s)=D(s)/C(s) ,
!       If IERROR=0: no errors detected
!                =1: invalid filter order l and k
!                                       in Chapter 7
!----------------------------------------------------------------------
    implicit none
    real*8 :: d(0:4),c(0:4)
    integer :: l,k,ln,ierror
    
    real*8 :: pi,orderk
    integer :: i
    
    pi = 4.0*atan(1.0)
    ierror = 1
    if ((l <= 0) .or. (k > int((l+1)/2))) return
    
    ierror = 0
    d = 0.0
    c = 0.0
    d(0) = 1.0
    c(0) = 1.0
    
    orderk = k-(l+1.0)/2.0
    if (orderk == 0.0) then
        c(1) = 1.0
    else
        c(1) = -2.0*cos((2*k+l-1)*pi/(2*l))
        c(2) = 1.0
    end if
    
    return
    end subroutine butwcf
    
