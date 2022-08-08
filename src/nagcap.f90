    program main
    
    use sacio
    use load_data
    use omp_lib
    
    ! inversion parameters
    character(len=40) :: filename,ddir,gdir
    character*80 :: weight
    integer :: leng(2),lenstf,tshift(2),pre_t(2)
    real*8 :: flh_p(2),flh_s(2),depth(3),lb(5),hb(5)
    real*8 :: dt,evdp
    real*8,allocatable :: stf(:)
    
    character*80 :: fnm,gdir0,outdir,newdir
    
    type(Seis),allocatable :: ob(:), g(:)
    real*8,allocatable :: rx(:),ry(:)
    character*6,allocatable :: stnm(:)
    integer :: nst
    
    real*8,allocatable :: vf(:,:),vm(:)
    
    real*8 :: err,sca_M,moment(6),mm(6),Mw
    integer :: idx(1),mexp
    real*8 :: eta,kai,iso,dc,clvd,pi
    
    interface
        function triangle(lenstf,shift)
            integer :: lenstf,shift
            real*8 :: triangle(shift+lenstf+1)
        end function triangle
        
        subroutine fitness(fault,ob,g,nst,dt,pre_t,leng,tshift,err,sca_M)
            use load_data
            real*8,intent(in) :: fault(5)
            integer,intent(in) :: nst,pre_t(2),leng(2),tshift(2)
            type(Seis),intent(in) :: ob(:),g(:)
            real*8,intent(in) :: dt
            real*8,intent(out) :: err,sca_M
        end subroutine fitness
        
        subroutine NA(ob,g,nst,dt,pre_t,leng,tshift,ns,nr,niter,nd,lb,hb,vm,vf)
            use load_data
            integer,intent(in) :: nst,ns,nr,niter,nd
            integer,intent(in) :: pre_t(2),leng(2),tshift(2)
            type(Seis),intent(in) :: ob(:),g(:)
            real*8,intent(in) :: dt
            real*8,intent(in) :: lb(nd),hb(nd)
            real*8,intent(out) :: vf(ns+ns*niter,nd),vm(ns+ns*niter)
        end subroutine NA
       
		subroutine deloc(ob,g,nst,rx,ry,stnm)
            use load_data
            integer,intent(in) :: nst
            type(Seis),allocatable :: ob(:),g(:)
            real*8,allocatable :: rx(:),ry(:)
            character*6,allocatable :: stnm(:)
        end subroutine deloc
 
    end interface
    
    pi = 4.0*atan(1.0)
    
    call getarg(1,newdir)    
    filename = 'input1.par'
!    call get_par('201906171455430/'//filename,ddir,gdir,leng,flh_p,flh_s,weight,lenstf,depth,tshift,lb,hb)
    call get_par(trim(newdir)//'/'//filename,ddir,gdir,leng,flh_p,flh_s,weight,lenstf,depth,tshift,lb,hb)

    outdir = gdir(1:scan(gdir,'/',.true.))
    
    allocate(stf(lenstf+0+1))
    stf = triangle(lenstf,0)
    
    lb = lb*pi/180.0
    hb = hb*pi/180.0
    ns = 1000
    nr = 1000
    nd = size(lb)
    niter = 50
    
    pre_t = int(leng*0.3)
    
    do evdp=depth(1),depth(2),depth(3)
        write(fnm,'(i)')int(evdp)
        write(*,*)trim(gdir)//adjustl(fnm)
        gdir0 = trim(trim(gdir)//adjustl(fnm))
        call loaddata(ddir,trim(gdir0),weight,leng,flh_p,flh_s,stf,evdp,ob,g,dt,rx,ry,stnm)
        
        allocate(vf(ns+ns*niter,nd))
        allocate(vm(ns+ns*niter))
        vf = 0.0
        vm = 0.0
        nst = size(ob)
        call NA(ob,g,nst,dt,pre_t,leng,tshift,ns,nr,niter,nd,lb,hb,vm,vf)
        
        open(10,file=gdir(1:scan(gdir,'/',.true.))//'result_'//trim(adjustl(fnm))//'.txt')
        do i=1,ns+ns*niter
            write(10,*)(vf(i,j),j=1,nd),vm(i)
        end do
        close(10)
        
        idx = minloc(vm)
        call fitness(vf(idx(1),:),ob,g,nst,dt,pre_t,leng,tshift,err,sca_M)
        Mw = log10(sca_M)/1.5-6.033
        
        call nmtensor(vf(idx(1),4),vf(idx(1),5),vf(idx(1),1),vf(idx(1),2),vf(idx(1),3),moment)
        mm(1)=moment(6);mm(2)=moment(1);mm(3)=moment(4);mm(4)=moment(3);mm(5)=-moment(5);mm(6)=-moment(2);

        eta = cos(vf(idx(1),4))
        kai = sin(vf(idx(1),5))
        iso = sign(1.0,eta)*eta*eta
        dc = (1.0-eta*eta)*(1-kai*kai)
        clvd = sign(1.0,kai)*(1.0-eta*eta)*kai*kai
        mexp = anint(log10(sca_M))

        write(fnm,'(i)')int(evdp)
        fnm = trim(outdir)//'output_'//trim(adjustl(fnm))//'.txt'
        print *,fnm
        open(10,file=fnm)
        write(10,*)gdir0
        write(10,'("Mw = ",f9.5)')Mw
        write(10,'("FM =",3(f9.4,1x),f8.4,1x,f8.4)')vf(idx(1),:)*180.0/pi
        write(10,'("M = ",6(f10.6,1x))')mm*(sca_M/(10.0**mexp))
        write(10,'("exp = ",i5)')mexp
        write(10,'("ISO = ",f8.4," DC = ",f8.4," CLVD = ",f8.4)')iso,dc,clvd
        write(10,'("ERR  = ",e11.4)')vm(idx(1))
        close(10)
        
        deallocate(vf)
        deallocate(vm)
		call deloc(ob,g,nst,rx,ry,stnm)
    end do
    
    deallocate(stf)
    
    stop
    end

    subroutine deloc(ob,g,nst,rx,ry,stnm)
    use load_data
    implicit none
    
    integer,intent(in) :: nst
    type(Seis),allocatable :: ob(:),g(:)
    real*8,allocatable :: rx(:),ry(:)
    character*6,allocatable :: stnm(:)
    
    integer :: ist
    
    deallocate(rx)
    deallocate(ry)
    deallocate(stnm)
    do ist=1,nst
        deallocate(ob(ist)%comp)
        deallocate(ob(ist)%pr)
        deallocate(ob(ist)%pz)
        deallocate(ob(ist)%sr)
        deallocate(ob(ist)%st)
        deallocate(ob(ist)%sz)
        
        deallocate(g(ist)%comp)
        deallocate(g(ist)%pr)
        deallocate(g(ist)%pz)
        deallocate(g(ist)%sr)
        deallocate(g(ist)%st)
        deallocate(g(ist)%sz)
    end do
    deallocate(ob)
    deallocate(g)
    return
    end subroutine deloc
