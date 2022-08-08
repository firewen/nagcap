    subroutine gen_p_data(fault,ob,g,nst,dt,pre_t,leng,tshift,sca_M,stnm,outdir)
    ! generate synthetic and observed data for plot
    use load_data
    implicit none
    
    real*8,intent(in) :: fault(5)
    integer,intent(in) :: nst,pre_t(2),leng(2),tshift(2)
    !type(Seis),intent(in) :: ob(nst),g(nst)
    type(Seis),intent(in) :: ob(:),g(:)
    real*8,intent(in) :: dt
    real*8,intent(in) :: sca_M
    character*6,allocatable :: stnm(:)
    character(len=*) :: outdir
    
    character*80 :: fnm
    
    real*8 :: strike,dip,rake,beta,gamma
    real*8,allocatable :: syn_pr(:,:),ob1_pr(:,:)
    real*8,allocatable :: syn_pz(:,:),ob1_pz(:,:)
    real*8,allocatable :: syn_sr(:,:),ob1_sr(:,:)
    real*8,allocatable :: syn_st(:,:),ob1_st(:,:)
    real*8,allocatable :: syn_sz(:,:),ob1_sz(:,:)
    
    real*8 :: cc_pr,cc_pz,cc_sr,cc_sz,cc_st
    integer :: is_pr,is_pz,is_sr,is_st,is_sz
    integer :: counter,ist,i
    
    real*8 :: moment(6)
    real*8 :: norm_ob,norm_syn
    real*8 :: corrcoef,misfit,norm
    real*8 :: maxamp,obmax,synmax
    
    strike = fault(1)
    dip = fault(2)
    rake = fault(3)
    beta = fault(4)
    gamma = fault(5)
    
    is_pr = 0;is_pz = 0;is_sr = 0;is_st = 0;is_sz = 0
    allocate(syn_pr(leng(1),nst));syn_pr = 0.0
    allocate(ob1_pr(leng(1),nst));ob1_pr = 0.0
    allocate(syn_pz(leng(1),nst));syn_pz = 0.0
    allocate(ob1_pz(leng(1),nst));ob1_pz = 0.0
    allocate(syn_sr(leng(2),nst));syn_sr = 0.0
    allocate(ob1_sr(leng(2),nst));ob1_sr = 0.0
    allocate(syn_st(leng(2),nst));syn_st = 0.0
    allocate(ob1_st(leng(2),nst));ob1_st = 0.0
    allocate(syn_sz(leng(2),nst));syn_sz = 0.0
    allocate(ob1_sz(leng(2),nst));ob1_sz = 0.0
    
    call nmtensor(beta,gamma,strike,dip,rake,moment)
    
    ! generate synthetic seismograms and align with observation, then compute the correlation coefficience
    counter = 0

    fnm=trim(outdir)//'/plotd.txt'
    open(10,file=fnm)
    do ist=1,nst
        print *, ist
        if(any(ob(ist).comp /=0)) then
            counter = counter+1
        else
            cycle
        end if
        
        norm_ob = 0.0
        norm_syn = 0.0
        
        ! generate one group of synthetic records
        if (ob(ist).comp(1) /= 0) then      ! Pr
            call align_seis(ob(ist),g(ist),nst,moment,'P','r',dt,pre_t(1),leng(1),tshift(1),    &
                            is_pr,ob1_pr(:,ist),syn_pr(:,ist))    
        end if
        
        if (ob(ist).comp(2) /= 0) then      ! Pz
            call align_seis(ob(ist),g(ist),nst,moment,'P','z',dt,pre_t(1),leng(1),tshift(1),    &
                            is_pz,ob1_pz(:,ist),syn_pz(:,ist))    
        end if
        
        if ((ob(ist).comp(1) /=0) .and. (ob(ist).comp(2) /= 0)) then
        ! if use Pr and Pz to invert, the time shift shoulb be same for Pr and Pz
            cc_pr = corrcoef(ob1_pr(:,ist),syn_pr(:,ist),leng(1))
            cc_pz = corrcoef(ob1_pz(:,ist),syn_pz(:,ist),leng(1))
            if (cc_pr > cc_pz) then
                call align_seis(ob(ist),g(ist),nst,moment,'P','z',dt,pre_t(1),leng(1),0,    &
                            is_pr,ob1_pz(:,ist),syn_pz(:,ist))
                is_pz = is_pr
            else
                call align_seis(ob(ist),g(ist),nst,moment,'P','r',dt,pre_t(1),leng(1),0,    &
                            is_pz,ob1_pr(:,ist),syn_pr(:,ist))
                is_pr = is_pz
            end if
        end if
        
        if (ob(ist).comp(3) /= 0) then      ! Sr
            call align_seis(ob(ist),g(ist),nst,moment,'S','r',dt,pre_t(2),leng(2),tshift(2),    &
                            is_sr,ob1_sr(:,ist),syn_sr(:,ist))    
        end if
        
        if (ob(ist).comp(4) /= 0) then      ! St
            call align_seis(ob(ist),g(ist),nst,moment,'S','t',dt,pre_t(2),leng(2),tshift(2),    &
                            is_st,ob1_st(:,ist),syn_st(:,ist))    
        end if
        
        if (ob(ist).comp(5) /= 0) then      ! Sz
            call align_seis(ob(ist),g(ist),nst,moment,'S','z',dt,pre_t(2),leng(2),tshift(2),    &
                            is_sz,ob1_sz(:,ist),syn_sz(:,ist))    
        end if
        
        if ((ob(ist).comp(3) /=0) .and. (ob(ist).comp(5) /= 0)) then
        ! if use Sr and Sz to invert, the time shift shoulb be same for Sr and Sz
            cc_sr = corrcoef(ob1_sr(:,ist),syn_sr(:,ist),leng(2))
            cc_sz = corrcoef(ob1_sz(:,ist),syn_sz(:,ist),leng(2))
            if (cc_sr > cc_sz) then
                call align_seis(ob(ist),g(ist),nst,moment,'S','z',dt,pre_t(2),leng(2),0,    &
                            is_sr,ob1_sz(:,ist),syn_sz(:,ist))
                is_sz = is_sr
            else
                call align_seis(ob(ist),g(ist),nst,moment,'S','r',dt,pre_t(2),leng(2),0,    &
                            is_sz,ob1_sr(:,ist),syn_sr(:,ist))
                is_sr = is_sz
            end if
        end if
        
        print *, is_pr,is_pz,is_sr,is_st,is_sz
        ! output epicenter distance, azimuth
        write(10,'(a6,1x,f8.4,1x,f8.4\)')stnm(ist), ob(ist).dist, ob(ist).az
        if (ob(ist).comp(1) /= 0) then ! Pr
            cc_pr = corrcoef(ob1_pr(:,ist),syn_pr(:,ist),leng(1))

            obmax = maxval(abs(ob1_pr(:,ist)))
            synmax = maxval(abs(syn_pr(:,ist)))*sca_M
            maxamp = max(obmax,synmax)
            ! output cc, time shift for pr
            write(10,'(1x,f8.4,1x,f8.4,1x,e8.1,1x,e8.1\)')cc_pr*100, is_pr*dt, obmax, synmax
            
            open(20,file=trim(outdir)//'/'//trim(stnm(ist))//'.pr')
            do i=1,leng(1)
                write(20,'(f8.4,1x,f8.4,1x,f8.4)')(i-1)*dt,ob1_pr(i,ist)/maxamp,syn_pr(i,ist)*sca_M/maxamp
            end do
            close(20)
        end if
        
        if (ob(ist).comp(2) /= 0) then ! Pz
            cc_pz = corrcoef(ob1_pz(:,ist),syn_pz(:,ist),leng(1))
            
            obmax = maxval(abs(ob1_pz(:,ist)))
            synmax = maxval(abs(syn_pz(:,ist)))*sca_M
            maxamp = max(obmax,synmax)
            ! output cc, time shift for pz
            write(10,'(1x,f8.4,1x,f8.4,1x,e8.1,1x,e8.1\)')cc_pz*100, is_pz*dt, obmax, synmax
            
            open(20,file=trim(outdir)//'/'//trim(stnm(ist))//'.pz')
            do i=1,leng(1)
                write(20,'(f8.4,1x,f8.4,1x,f8.4)')(i-1)*dt,ob1_pz(i,ist)/maxamp,syn_pz(i,ist)*sca_M/maxamp
            end do
            close(20)
        end if
        
        if (ob(ist).comp(3) /= 0) then ! Sr
            cc_sr = corrcoef(ob1_sr(:,ist),syn_sr(:,ist),leng(2))
            
            obmax = maxval(abs(ob1_sr(:,ist)))
            synmax = maxval(abs(syn_sr(:,ist)))*sca_M
            maxamp = max(obmax,synmax)
            ! output cc, time shift for sr
            write(10,'(1x,f8.4,1x,f8.4,1x,e8.1,1x,e8.1\)')cc_sr*100, is_sr*dt, obmax, synmax
            
            open(20,file=trim(outdir)//'/'//trim(stnm(ist))//'.sr')
            do i=1,leng(2)
                write(20,'(f8.4,1x,f8.4,1x,f8.4)')(i-1)*dt,ob1_sr(i,ist)/maxamp,syn_sr(i,ist)*sca_M/maxamp
            end do
            close(20)
        end if
        
        if (ob(ist).comp(4) /= 0) then ! St
            cc_st = corrcoef(ob1_st(:,ist),syn_st(:,ist),leng(2))
            
            obmax = maxval(abs(ob1_st(:,ist)))
            synmax = maxval(abs(syn_st(:,ist)))*sca_M
            maxamp = max(obmax,synmax)
            ! output cc, time shift for st
            write(10,'(1x,f8.4,1x,f8.4,1x,e8.1,1x,e8.1\)')cc_st*100, is_st*dt, obmax, synmax
            
            open(20,file=trim(outdir)//'/'//trim(stnm(ist))//'.st')
            do i=1,leng(2)
                write(20,'(f8.4,1x,f8.4,1x,f8.4)')(i-1)*dt,ob1_st(i,ist)/maxamp,syn_st(i,ist)*sca_M/maxamp
            end do
            close(20)
        end if
        
        if (ob(ist).comp(5) /= 0) then ! Sz
            cc_sz = corrcoef(ob1_sz(:,ist),syn_sz(:,ist),leng(2))
            
            obmax = maxval(abs(ob1_sz(:,ist)))
            synmax = maxval(abs(syn_sz(:,ist)))*sca_M
            maxamp = max(obmax,synmax)
            ! output cc, time shift for sz
            write(10,'(1x,f8.4,1x,f8.4,1x,e8.1,1x,e8.1\)')cc_sz*100, is_sz*dt, obmax, synmax
            
            open(20,file=trim(outdir)//'/'//trim(stnm(ist))//'.sz')
            do i=1,leng(2)
                write(20,'(f8.4,1x,f8.4,1x,f8.4)')(i-1)*dt,ob1_sz(i,ist)/maxamp,syn_sz(i,ist)*sca_M/maxamp
            end do
            close(20)
        end if
        write(10,'(1x)')
        
        
    end do    
    close(10)
    
    deallocate(syn_pr)
    deallocate(ob1_pr)
    deallocate(syn_pz)
    deallocate(ob1_pz)
    deallocate(syn_sr)
    deallocate(ob1_sr)
    deallocate(syn_st)
    deallocate(ob1_st)
    deallocate(syn_sz)
    deallocate(ob1_sz)
    return
    
    end subroutine gen_p_data