    subroutine align_seis(ob,g,nst,moment,wtype,cmp,dt,pre_t,leng,tshift,ishift,ob1,syn)
    use load_data
    implicit none
    
    integer,intent(in) :: nst,pre_t,leng,tshift
    type(Seis),intent(in) :: ob,g
    real*8,intent(in) :: dt
    real*8,intent(in) :: moment(6)
    character,intent(in) :: wtype,cmp
    integer,intent(inout) :: ishift
    real*8,intent(out) :: ob1(leng),syn(leng)
    
    real*8 :: T,B
    real*8 :: g0(leng,6)
    real*8,allocatable :: ob0(:),ob2(:)
    integer :: nt,cut_start,cut_end
    real*8 :: tmp,oldcc,tmpob(leng)
    real*8 :: acor(2*leng-1)
    integer :: lag(2*leng-1)
    integer :: idx
    integer :: idx0(1)
    integer :: i
    real*8 :: corrcoef
    
    if ((wtype == 'P') .and. (cmp == 'r')) then
        T = g.T1
        B = ob.rb
        nt = size(ob.pr)
        allocate(ob0(nt))
        ob0 = reshape(ob.pr,(/nt/))
        g0 = g.pr
    end if
    if ((wtype == 'P') .and. (cmp == 'z')) then
        T = g.T1
        B = ob.zb
        nt = size(ob.pz)
        allocate(ob0(nt))
        ob0 = reshape(ob.pz,(/nt/))
        g0 = g.pz
    end if
    if ((wtype == 'S') .and. (cmp == 'r')) then
        T = g.T2
        B = ob.rb
        nt = size(ob.sr)
        allocate(ob0(nt))
        ob0 = reshape(ob.sr,(/nt/))
        g0 = g.sr
    end if
    if ((wtype == 'S') .and. (cmp == 't')) then
        T = g.T2
        B = ob.tb
        nt = size(ob.st)
        allocate(ob0(nt))
        ob0 = reshape(ob.st,(/nt/))
        g0 = g.st
    end if
    if ((wtype == 'S') .and. (cmp == 'z')) then
        T = g.T2
        B = ob.zb
        nt = size(ob.sz)
        allocate(ob0(nt))
        ob0 = reshape(ob.sz,(/nt/))
        g0 = g.sz
    end if
    
    if (tshift == 0) then
        if (ishift == 0) then
        ! there is no time shift, directly cut data according to phase arrivals
            cut_start = int((T-B)/dt)-pre_t
            cut_end = cut_start+leng-1
        ! referred to theoretical phase arrivals to time shift data
        ! in order to keep time shift for Pr and Pz or time shift for Sr, Sz and St the same 
        else if (ishift > 0) then
            cut_start = int((T-B)/dt)-pre_t
            cut_start = cut_start-ishift
            cut_end = cut_start+leng-1
        else if (ishift < 0) then
            cut_start = int((T-B)/dt)-pre_t
            cut_start = cut_start-ishift
            cut_end = cut_start+leng-1
        end if
        
        if (cut_start < 0) then
            cut_start = 1
            cut_end = leng
        end if
        if (cut_end > nt) then
            allocate(ob2(cut_end))
            ob2 = 0.0
            ob2(1:nt) = ob0
        else
            allocate(ob2(nt))
            ob2 = ob0
        end if
        
        ob1 = ob2(cut_start:cut_end)
        syn = reshape(matmul(g0,reshape(moment,(/6,1/))),(/leng/))
        
        call taper(ob1,leng,pre_t)
        call taper(syn,leng,pre_t)
        
        deallocate(ob0)
        deallocate(ob2)
        return
    else
        cut_start = int((T-B)/dt)-pre_t
        cut_end = cut_start+leng-1
        
        if (cut_start < 0) then
            cut_start = 1
            cut_end = leng
        end if
        ob1 = ob0(cut_start:cut_end)
        
        syn = reshape(matmul(g0,reshape(moment,(/6,1/))),(/leng/))
        
        tmp = corrcoef(syn,ob1,leng)
        oldcc = tmp
        
        call xcorrf(syn,ob1,leng,acor,lag) 
        !print *, 'idx0 = ',maxloc(acor)
        call findpeaks(acor,lag,2*leng-1,tshift,idx)
        !idx0 = maxloc(acor)
        !idx = idx0(1)
        
        if (lag(idx) < 0) then
            if ((lag(idx) <= -tshift) .or. (lag(idx) >= 0)) then
                ishift = 0
                tmpob = ob1
            else
                cut_start = cut_start-lag(idx)
                cut_end = cut_start+leng-1
                
                if (cut_end > nt) then
                    tmpob = 0.0
                    tmpob(1:nt-cut_start+1) = ob0(cut_start:nt)
                else
                    tmpob = ob0(cut_start:cut_end)
                end if
                ishift = lag(idx)
                
                tmp = corrcoef(syn,tmpob,leng)
                if (tmp < oldcc) then
                    ishift = 0
                    tmpob = ob1
                end if
            end if
        else
            if ((lag(idx) <= 0) .or. (lag(idx) >= tshift)) then
                ishift = 0
                tmpob = ob1
            else
                cut_start = cut_start-lag(idx)
                cut_end = cut_start+leng-1  
                tmpob = ob0(cut_start:cut_end)
                ishift = lag(idx)
                
                tmp = corrcoef(syn,tmpob,leng)
                if (tmp < oldcc) then
                    ishift = 0
                    tmpob = ob1
                end if
            end if
        end if
        
        ob1 = tmpob

        call taper(ob1,leng,pre_t)
        call taper(syn,leng,pre_t)
        
        deallocate(ob0)
        return
    end if
    
    end subroutine align_seis
    
    subroutine taper(x,nt,pre_t)
    implicit none
    
    integer,intent(in) :: nt,pre_t
    real*8,intent(inout) :: x(nt)
    
    real*8 :: pi,omega,f0,f1
    integer :: i
    
    pi = 4.0*atan(1.0)
    omega = pi/pre_t
    f0 = 0.5
    f1 = 0.5
    
    do i=1,pre_t
        x(i) = x(i)*(f0-f1*cos(omega*(i-1)))
    end do
    
    return
    end subroutine taper
    
    function corrcoef(x,y,n)
    implicit none
    
    integer,intent(in) :: n
    real*8,intent(in) :: x(n,1),y(n,1)
    real*8 :: corrcoef
    
    real*8 :: xmean,ymean,x2sum(1,1),y2sum(1,1),xysum(1,1)
    
    xmean = sum(x)/n
    ymean = sum(y)/n
    
    x2sum = matmul(transpose(x-xmean),(x-xmean))
    y2sum = matmul(transpose(y-ymean),(y-ymean))
    xysum = matmul(transpose(x-xmean),(y-ymean))
    
    corrcoef = xysum(1,1)/sqrt(x2sum(1,1)*y2sum(1,1))
    
    return
    end function corrcoef
    
    subroutine xcorrf(x,y,n,z,lag)
    implicit none
    
    integer,intent(in) :: n
    real*8,intent(in) :: x(n),y(n)
    real*8,intent(out) :: z(2*n-1)
    integer,intent(out) :: lag(2*n-1)
    
    integer :: i,m,leng
    real*8,allocatable :: xx(:),yy(:)
    complex*16,allocatable :: xspec(:),yspec(:),zspec(:)
    
    logical :: ispower_of_two
    integer :: numTo2N
    
    do i=1-n,n-1
        lag(i+n) = i
    end do
    
    if (.not. ispower_of_two(2*n-1)) then
        leng = numTo2N(2*n-1)
        m = log(real(leng))/log(2.0)
    end if
    
    allocate(xx(leng))
    allocate(yy(leng))
    xx = 0.0
    yy = 0.0
    xx(1:n) = x
    yy(1:n) = y
    
    allocate(xspec(leng))
    allocate(yspec(leng))
    allocate(zspec(leng))
    xspec = cmplx(xx,0)
    yspec = cmplx(yy,0)
    deallocate(xx)
    deallocate(yy)
    
    call fft(xspec,m,1)
    call fft(yspec,m,1)
    
    zspec = conjg(yspec)*xspec
    
    call fft(zspec,m,-1)

    zspec = cshift(real(zspec),1-n)
    z = real(zspec(1:2*n-1))

    !deallocate(xx)
    !deallocate(yy)
    deallocate(xspec)
    deallocate(yspec)
    deallocate(zspec)
    return
    end subroutine xcorrf
    
    function ispower_of_two(a) result(b)
    implicit none
    integer,intent(in) :: a
    logical :: b
    if((iand(a,-a) == a) .and. (a /= 0)) then
        b = .true.
    else
        b = .false.
    end if
    end function ispower_of_two
    
    function numTo2N(a) result(i)
    implicit none
    integer,intent(in) :: a
    
    integer :: i,tmp
    
    i = 1
    tmp = a
    
    do 
        tmp = int(tmp/2)
        i = shiftl(i,1)
        if (tmp == 0) exit
    end do
    
    if (i < a) then
        i = shiftl(i,1)
    end if
    
    return
    end function numTo2N        

    subroutine fft(a,m,flag)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      按时间抽取的FFT子程序:
!c      a: 复数数组，长度为2**m，存放要进行fft的序列；
!c      m: 序列长度为 2**m.
!c      flag: FFT或IFFT的选择参数：
!c            flag=1  -------  FFT;
!c            flag=-1 ------- IFFT.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

    implicit none
    integer    l,m,n,n2,i,j,k,l1,l2,ip,flag
    real*8       pi
    complex*16   a(2**m),u,w,t

    pi=4d0*atan(1d0)
    n=2**m                      !   抽样点数
    n2=n/2            

    if(flag.eq.-1) then         ! 
        do i=1,n                !
            a(i)=conjg(a(i))    !   做ifft时对数组a取共轭
        end do               
    end if                 

!ccccccccccccccccccc  倒位序程序段  ccccccccccccccccccccc

    j=1
    do i=1,n-1
        if(i.lt.j) then     ! 
            t=a(j)          !
            a(j)=a(i)       !  如果i<j，互换a(i)和a(i)
            a(i)=t          !
        end if              !

        k=n2                !     
        do while(k.lt.j)    !      
            j=j-k           !  已知j,用反向加法
            k=k/2           !  求下一个倒位序号
        end do              !  
        j=j+k               !
    end do
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    do l=1,m                 ! 第一层循环，l为级数，从1到m
        l1=2**l              ! l1 : 第l级中的乘数个数 
        l2=l1/2              ! l2 : 第l级的蝶形间隔
        w=cmplx(cos(pi/float(l2)),-sin(pi/float(l2)))
	                         ! w  : 第l级的乘数底数
        u=cmplx(1.0,0.0)
        do j=1,l2            ! 第二层循环，乘数控制
            do i=j,n,l1      ! 第三层循环，群控制
                ip=i+l2      ! i:蝶形左上角序号，ip:左下脚序号
                t=a(ip)*u
                a(ip)=a(i)-t ! 蝶形运算
                a(i)=a(i)+t
            end do
            u=u*w            ! 第二层循环中调整乘数
        end do 
    end do

    if(flag.eq.-1) then      ! 做ifft时对数组再a取共轭，并除以n
        do i=1,n
            a(i)=conjg(a(i))/n
        end do
    end if     

    return
    end subroutine fft