module load_data
    public :: loaddata
    
    type :: Seis
        real*8,allocatable :: pr(:,:),pz(:,:),sr(:,:),st(:,:),sz(:,:),comp(:)
        real*8 :: dist,az,rb,tb,zb,t1,t2
    end type Seis
    
    !integer :: nst,leng(2)
    !type(Seis),allocatable :: ob(:), g(:)
    
    contains
    subroutine loaddata(ddir,gdir,weight,leng,flh_p,flh_s,stf,evdp,ob,g,dt,rx,ry,stnm)
    ! Input
    ! ddir : the directory for observation data
    ! gdir : the directory for empirical green functions
    ! weight : the weight file for each component of different station
    ! leng : the length of data used in inversion leng(1) for P wave, 
    !       leng(2) for S wave
    ! flh_p : the parameters for a bandpass filter for P wave, [lf hf];
    ! flh_s : the parameters for a bandpass filter for S wave, [lf hf];
    ! stf : the source time function (default is a triangle);
    ! Output:
    ! ob : 3 components data
    ! g  : empirical green functions for each station with 3 components
    ! dt : sampling time
    use sacio
    use raytracing
    implicit none
    
    integer :: nst,leng(2)
    real*8,allocatable :: rx(:),ry(:),rz(:),cmp(:,:),stf(:)
    character*6,allocatable :: stnm(:)
    character(len=*) :: ddir,gdir,weight
    real*8 :: flh_p(2),flh_s(2),evdp
    
    real*8 :: pi,sx,sy,r0,fai
    
    integer :: i,j,ist,ierror
    integer :: pre_t(2), nt
    real*8 :: dt, T1, T2, rB, tB, zB, dist, az, dto
    integer :: status

    ! the parameter for butterworth filter, a order 3 filter has been used
    real*8,allocatable :: fb_p(:),fa_p(:),fb_s(:),fa_s(:)
    real*8,allocatable :: frq(:)    
    real*8 :: srate     ! sampling rate
    integer :: lp,ls,iband
    
    type(Seis),allocatable :: ob(:), g(:)
    real*8,allocatable :: obr(:),obt(:),obz(:),tmp(:)
    character(len=80) :: buffer
    
    type(sachead) :: head
    real,allocatable :: dd(:)
    real*8 :: dtg
    real*8,allocatable :: ifun(:,:),tsac(:)
    
    ! velocity model
    type(Velocity) :: vmodel
    integer :: nlayer
    character*40 :: vfnm
    real*8 :: rho,fp
        
    interface
        subroutine butter(frq,fs,iband,l,bn,an)
            real*8,allocatable :: frq(:)
            real*8 :: fs
            integer :: l,iband
            real*8,allocatable :: bn(:),an(:)
        end subroutine butter
    end interface
    
    pi = 4.0*atan(1.0)
    
    ! load weight file
    open(10,file=weight)
    nst = 0
    do while(.true.)
        read(10,fmt='(a79)',iostat=status) buffer
        if(status /= 0) exit
        nst = nst+1
    end do
    close(10)

    allocate(stnm(nst))
    allocate(rx(nst))
    allocate(ry(nst))
    allocate(rz(nst))
    allocate(cmp(nst,5))
    open(10,file=weight)
    do ist=1,nst
        read(10,*) stnm(ist),rx(ist),ry(ist),rz(ist),(cmp(ist,j),j=1,5)
    end do
    close(10)
    
    ! load velocity model file
    vfnm = 'layered_model.dat'
    call model_info(weight(1:scan(weight,'/',.true.))//vfnm,vmodel)
    
    sx = 0.0; sy = 0.0
    
    allocate(ob(nst))
    allocate(g(nst))
    
    pre_t = int(leng*0.3)    ! pre_t(1) for P; pre_t(2) for S

    do ist=1,nst
        ! observation
        if ((cmp(ist,1) /= 0) .or. (cmp(ist,3) /= 0)) then
            call getsac(trim(ddir)//stnm(ist)//'.r',obr,nt,dt,T1,T2,rB,dist,az)
        end if

        if (cmp(ist,4) /= 0) then
            call getsac(trim(ddir)//stnm(ist)//'.t',obt,nt,dt,T1,T2,tB,dist,az)
        end if

        if ((cmp(ist,2) /= 0) .or. (cmp(ist,5) /= 0)) then
            call getsac(trim(ddir)//stnm(ist)//'.z',obz,nt,dt,T1,T2,zB,dist,az)
        end if

        if (nt < leng(1)) then
            print *, 'the length of the used data should be larger than the length of the observed data'
            return
        end if

        dto = dt

        allocate(ob(ist)%comp(5))
        allocate(ob(ist)%pr(nt,1))
        allocate(ob(ist)%pz(nt,1))
        allocate(ob(ist)%sr(nt,1))
        allocate(ob(ist)%st(nt,1))
        allocate(ob(ist)%sz(nt,1))

        ob(ist)%dist = dist
        ob(ist)%az = az
        srate = 1.0/dt
        
        ! generate filter parameters
        lp = 3
        if (flh_p(1) == 0) then
            iband = 1   !lp
            allocate(frq(1))
            frq(1) = flh_p(2)
            call butter(frq,srate,iband,lp,fb_p,fa_p)
        else
            iband = 3   ! bp
            allocate(frq(2))
            frq(1) = flh_p(1)
            frq(2) = flh_p(2)
            call butter(frq,srate,iband,lp,fb_p,fa_p)
        end if
        deallocate(frq)
        ls =3
        if (flh_s(1) == 0) then
            iband = 1   !lp
            allocate(frq(1))
            frq(1) = flh_s(2)
            call butter(frq,srate,iband,ls,fb_s,fa_s)
        else
            iband = 3   ! bp
            allocate(frq(2))
            frq(1) = flh_s(1)
            frq(2) = flh_s(2)
            call butter(frq,srate,iband,ls,fb_s,fa_s)
        end if
        deallocate(frq)
        
        ! for Pr
        if (cmp(ist,1) /= 0) then
            call fitout(fb_p,fa_p(1:lp),lp,lp,obr,nt,ob(ist).pr,ierror)
            ob(ist).comp(1) = cmp(ist,1)
            ob(ist).rb = rB
        end if
        
        ! for Sr
        if (cmp(ist,3) /= 0) then
            call fitout(fb_s,fa_s(1:ls),ls,ls,obr,nt,ob(ist).sr,ierror)
            ob(ist).comp(3) = cmp(ist,3)
            ob(ist).rb = rB
        end if
        
        ! for St
        if (cmp(ist,4) /= 0) then
            call fitout(fb_s,fa_s(1:ls),ls,ls,obt,nt,ob(ist).st,ierror)
            ob(ist).comp(4) = cmp(ist,4)
            ob(ist).tb = tB
        end if
        
        ! for Pz
        if (cmp(ist,2) /= 0) then
            call fitout(fb_p,fa_p(1:lp),lp,lp,obz,nt,ob(ist).pz,ierror)
            ob(ist).comp(2) = cmp(ist,2)
            ob(ist).zb = zB
        end if
        
        ! for Sz
        if (cmp(ist,5) /= 0) then
            call fitout(fb_s,fa_s(1:ls),ls,ls,obz,nt,ob(ist).sz,ierror)
            ob(ist).comp(5) = cmp(ist,5)
            ob(ist).zb = zB
        end if
        
        deallocate(obr)
        deallocate(obt)
        deallocate(obz)
        
        ! Green functions
        allocate(g(ist)%comp(5))
        allocate(g(ist)%pr(leng(1),6))
        allocate(g(ist)%pz(leng(1),6))
        allocate(g(ist)%sr(leng(2),6))
        allocate(g(ist)%st(leng(2),6))
        allocate(g(ist)%sz(leng(2),6))
        
        write(buffer,'("g_"a6,i2.2".sac")')stnm(ist),1
        call sacio_readhead(trim(gdir)//'/'//trim(buffer), head, ierror)
        nt = head%npts
        
        dtg = head%delta
        if (abs(dto-dtg) > 1e-5) then
            write(*,*)'The dt of observation (',dto,') is not equal to the green functions(',dtg,')'
            return
        end if
        
        allocate(ifun(nt,10))
        ! load I component for green function and convolution with STF
        if ((cmp(ist,1) /= 0) .or. (cmp(ist,3) /= 0)) then      ! R
            do i=3,6
                write(buffer,'("g_"a6,i2.2".sac")')stnm(ist),i
                call sacio_readsac(trim(gdir)//'/'//trim(buffer), head, dd, ierror)
                
                allocate(tsac(nt))
                
                tsac = real(dd,kind=8)
                
                allocate(tmp(nt+size(stf)-1))
                
                call conv(tsac,stf,tmp,nt,size(stf),nt+size(stf)-1,ierror)
                ifun(:,i) = tmp(1:nt)
                
                deallocate(dd)
                deallocate(tsac)
                deallocate(tmp)
            end do
        end if
        
        if (cmp(ist,4) /= 0) then       ! T
            do i=1,2
                write(buffer,'("g_"a6,i2.2".sac")')stnm(ist),i
                call sacio_readsac(trim(gdir)//'/'//trim(buffer), head, dd, ierror)
                
                allocate(tsac(nt))
                
                tsac = real(dd,kind=8)
                
                allocate(tmp(nt+size(stf)-1))
                
                call conv(tsac,stf,tmp,nt,size(stf),nt+size(stf)-1,ierror)
                ifun(:,i) = tmp(1:nt)
                
                deallocate(dd)
                deallocate(tsac)
                deallocate(tmp)
            end do
        end if
        
        if ((cmp(ist,2) /= 0) .or. (cmp(ist,5) /= 0)) then      ! Z
            do i=7,10
                write(buffer,'("g_"a6,i2.2".sac")')stnm(ist),i
                call sacio_readsac(trim(gdir)//'/'//trim(buffer), head, dd, ierror)
                
                allocate(tsac(nt))
                
                tsac = real(dd,kind=8)
                
                allocate(tmp(nt+size(stf)-1))
                
                call conv(tsac,stf,tmp,nt,size(stf),nt+size(stf)-1,ierror)
                ifun(:,i) = tmp(1:nt)
                
                deallocate(dd)
                deallocate(tsac)
                deallocate(tmp)
            end do
        end if
        
        r0 = sqrt((rx(ist)-sx)**2+(ry(ist)-sy)**2)
        if (rx(ist) >= sx) then
            fai = acos((ry(ist)-sy)/r0)
        else
            fai = 2*pi-acos((ry(ist)-sy)/r0)
        end if
        
        call ray_out(vmodel,r0,evdp,'P',fp,t1)
        call ray_out(vmodel,r0,evdp,'S',fp,t2)
        
        ! for Pr green
        if (cmp(ist,1) /= 0) then
            call generate_g(t1,dt,pre_t(1),leng(1),fb_p,fa_p,lp,ifun,nt,fai,'R',g(ist).pr)
            g(ist).comp(1) = cmp(ist,1)
            g(ist).t1 = t1
        end if
        
        ! for Sr green
        if (cmp(ist,3) /= 0) then
            call generate_g(t2,dt,pre_t(2),leng(2),fb_s,fa_s,ls,ifun,nt,fai,'R',g(ist).sr)
            g(ist).comp(3) = cmp(ist,3)
            g(ist).t2 = t2
        end if
        
        ! for St green
        if (cmp(ist,4) /= 0) then
            call generate_g(t2,dt,pre_t(2),leng(2),fb_s,fa_s,ls,ifun,nt,fai,'T',g(ist).st)
            g(ist).comp(4) = cmp(ist,4)
            g(ist).t2 = t2
        end if
        
        ! the observation coordinate is ENU, so the upward vibration is positive.
        ! however the coordinate used in green function computation is NED and the
        ! downward vibration is minus.
        
        ! for Pz green
        if (cmp(ist,2) /= 0) then
            call generate_g(t1,dt,pre_t(1),leng(1),fb_p,fa_p,lp,ifun,nt,fai,'Z',g(ist).pz)
            g(ist).pz = -1.0*g(ist).pz
            g(ist).comp(2) = cmp(ist,2)
            g(ist).t1 = t1
        end if
        
        if (cmp(ist,5) /= 0) then
            call generate_g(t2,dt,pre_t(2),leng(2),fb_s,fa_s,ls,ifun,nt,fai,'Z',g(ist).sz)
            g(ist).sz = -1.0*g(ist).sz
            g(ist).comp(5) = cmp(ist,5)
            g(ist).t2 = t2
        end if
        
        deallocate(fb_p)
        deallocate(fa_p)
        deallocate(fb_s)
        deallocate(fa_s)
        deallocate(ifun)
        
    end do

    !deallocate(rx)
    !deallocate(ry)
    deallocate(rz)

    return
    end subroutine loaddata

    subroutine getsac(kname,x,npts,dt,T1,T2,B,dist,az)
    use sacio
    
    character(len=*) :: kname
    real*8,allocatable :: x(:)
    real*8 :: dt,T1,T2,B,dist,az
    integer :: npts
    
    integer :: flag
    type(sachead) :: head
    real,allocatable :: dd(:)
    
    interface
        function detrend(y,n)
            integer :: n
            real*8 :: y(n)
            real*8 :: detrend(n)
        end function detrend
    end interface
    
    call sacio_readhead(kname, head, flag)
    npts = head%npts
    allocate(x(npts),stat=flag)
    if (flag /= 0) then
        write(*, *) "Not enough memory for ifun"
        return
    end if
    
    !allocate(dd(npts),stat=flag)
    call sacio_readsac(kname, head, dd, flag)
    x = real(dd, kind=8)
    x = detrend(x,npts)
    x = x-sum(x)/real(npts,kind=8)
    dt = head%delta
    T1 = head%T1
    T2 = head%T2
    B = head%B
    dist = head%dist
    az = head%az

    deallocate(dd)
    
    return
    end subroutine getsac
    
    subroutine generate_g(t,dt,pre_t,leng,fb,fa,lo,ifun,nt,fai,dtype,green)
    implicit none
    integer :: pre_t,leng,lo,nt
    real*8 :: t,dt,fb(0:lo),fa(0:lo),ifun(nt,10),fai
    character :: dtype
    
    real*8 :: green(leng,6)
    
    real*8 :: ifun1(nt,10),ifun2(leng,10)
    real*8,allocatable :: ifun3(:,:)
    
    integer :: cut_start,cut_end
    real*8 :: B
    integer :: i,ierror
    
    !integer :: j
    
    B = 0.0
    cut_start = int((t-B)/dt)-pre_t
    cut_end = cut_start+leng-1
    
    do i=1,10
        call fitout(fb,fa(1:lo),lo,lo,ifun(:,i),nt,ifun1(:,i),ierror)
    end do
    
    allocate(ifun3(nt,10))
    ifun3 = ifun1
    if (cut_start <=  0) then
        cut_start = 1; cut_end = leng
    end if
    if (cut_end > nt) then
        deallocate(ifun3)
        allocate(ifun3(cut_end,10))
        ifun3 = 0.0
        ifun3(1:nt,:) = ifun1
    end if
    
    ifun2 = ifun3(cut_start:cut_end,:)
    
    !print *, shape(ifun1),shape(ifun2)
    
    call assemble_g(ifun2,leng,fai,dtype,green)

    !open(10,file='ifun.txt')
    !do i=1,leng
    !    write(10,'(10(e20.12,1x))') (ifun2(i,j),j=1,10)
    !end do
    !close(10)
    !open(10,file='green.txt')
    !do i=1,leng
    !    write(10,'(6(e20.12,1x))') (green(i,j),j=1,6)
    !end do
    !close(10)
    
    deallocate(ifun3)
    return
    end subroutine generate_g
    
    subroutine assemble_g(ifun,leng,fai,dtype,green)
    implicit none
    integer,intent(in) :: leng
    real*8,intent(in) :: ifun(leng,10),fai
    character,intent(in) :: dtype
    
    real*8,intent(out) :: green(leng,6)
    
    ! R
    if (dtype == 'R') then
        green(:,1) = ifun(:,3)*cos(2.0*fai)*0.5+ifun(:,5)*0.5
        green(:,2) = ifun(:,3)*sin(2.0*fai)
        green(:,3) = ifun(:,4)*cos(fai)
        green(:,4) =-ifun(:,3)*cos(2*fai)*0.5+ifun(:,5)*0.5
        green(:,5) = ifun(:,4)*sin(fai)
        green(:,6) = ifun(:,6)
    end if
    ! T
    if (dtype == 'T') then
        green(:,1) =-ifun(:,2)*sin(2.0*fai)*0.5
        green(:,2) = ifun(:,2)*cos(2.0*fai)
        green(:,3) = ifun(:,1)*sin(fai)
        green(:,4) = ifun(:,2)*sin(2.0*fai)*0.5
        green(:,5) =-ifun(:,1)*cos(fai) 
        green(:,6) = 0.0
    end if
    ! Z
    if (dtype == 'Z') then
        green(:,1) = ifun(:,7)*cos(2.0*fai)*0.5+ifun(:,9)*0.5
        green(:,2) = ifun(:,7)*sin(2.0*fai)
        green(:,3) = ifun(:,8)*cos(fai)
        green(:,4) =-ifun(:,7)*cos(2.0*fai)*0.5+ifun(:,9)*0.5
        green(:,5) = ifun(:,8)*sin(fai)
        green(:,6) = ifun(:,10)
    end if
    
    return
    end subroutine assemble_g
end module load_data
    

