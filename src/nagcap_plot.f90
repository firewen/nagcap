    program main
    
    use sacio
    use load_data
!    use omp_lib
    
    ! inversion parameters
    character(len=40) :: filename,ddir,gdir
    character*80 :: weight
    integer :: leng(2),lenstf,tshift(2),pre_t(2)
    real*8 :: flh_p(2),flh_s(2),depth(3),lb(5),hb(5)
    real*8 :: dt,evdp
    real*8,allocatable :: stf(:)
    
    character*80 :: fnm,gdir0,outdir,newdir,list
    
    type(Seis),allocatable :: ob(:), g(:)
    real*8,allocatable :: rx(:),ry(:)
    character*6,allocatable :: stnm(:)
    integer :: nst

    real*8,allocatable :: vf(:,:)
    real*8 :: sca_M,Mw
    
    interface
        function triangle(lenstf,shift)
            integer :: lenstf,shift
            real*8 :: triangle(shift+lenstf+1)
        end function triangle
        
        subroutine gen_p_data(fault,ob,g,nst,dt,pre_t,leng,tshift,sca_M,stnm,outdir)
            use load_data
            real*8,intent(in) :: fault(5)
            integer,intent(in) :: nst,pre_t(2),leng(2),tshift(2)
            type(Seis),intent(in) :: ob(:),g(:)
            real*8,intent(in) :: dt
            real*8,intent(in) :: sca_M
            character*6,allocatable :: stnm(:)
            character(len=*) :: outdir
        end subroutine gen_p_data
        
	subroutine deloc(ob,g,nst,rx,ry,stnm)
            use load_data
            integer,intent(in) :: nst
            type(Seis),allocatable :: ob(:),g(:)
            real*8,allocatable :: rx(:),ry(:)
            character*6,allocatable :: stnm(:)
        end subroutine deloc
 
    end interface
    
    pi = 4.0*atan(1.0)
    
!    list="10   "
    call getarg(2,list)    
    read(list,*)evdp
    print *, evdp
    call getarg(1,newdir)    
!    newdir='201906171455430'
    !newdir='./original'
    filename = 'input1.par'    
    
    call get_par(trim(newdir)//'/'//filename,ddir,gdir,leng,flh_p,flh_s,weight,lenstf,depth,tshift,lb,hb)
    nd = size(lb)
    
    outdir = gdir(1:scan(gdir,'/',.true.))
    
    allocate(stf(lenstf+0+1))
    stf = triangle(lenstf,0)
    
    pre_t = int(leng*0.3)
    !do evdp=depth(1),depth(2),depth(3)
    !evdp = 1
    !evdp = 10
        write(fnm,'(i)')int(evdp)
        write(*,*)trim(gdir)//adjustl(fnm)
        gdir0 = trim(trim(gdir)//adjustl(fnm))
        call loaddata(ddir,trim(gdir0),weight,leng,flh_p,flh_s,stf,evdp,ob,g,dt,rx,ry,stnm)
        nst = size(ob)
        
        allocate(vf(1,nd))
        fnm = trim(outdir)//'output_'//trim(adjustl(fnm))//'.txt'
        print *,fnm
        open(10,file=fnm)
        read(10,'(a80)') list
        read(10,'(a80)') list
        read(10,'("Mw = ",a80)')list
        read(list,*)Mw
        read(10,'("FM =",a80)')list
        read(list,*)vf(1,:)
        close(10)
        
        sca_M = 10**((Mw+6.033)*1.5)
        call gen_p_data(vf(1,:)*pi/180.0,ob,g,nst,dt,pre_t,leng,tshift,sca_M,stnm,trim(outdir)//'figd/')
        
        deallocate(vf)
	call deloc(ob,g,nst,rx,ry,stnm)
    !end do
    
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
