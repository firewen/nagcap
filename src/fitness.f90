    subroutine fitness(fault,ob,g,nst,dt,pre_t,leng,tshift,err,sca_M)
    !subroutine fitness(fault,err,sca_M)
    use load_data
    implicit none
    
    real*8,intent(in) :: fault(5)
    integer,intent(in) :: nst,pre_t(2),leng(2),tshift(2)
    !type(Seis),intent(in) :: ob(nst),g(nst)
    type(Seis),intent(in) :: ob(:),g(:)
    real*8,intent(in) :: dt
    real*8,intent(out) :: err,sca_M
    
    real*8 :: strike,dip,rake,beta,gamma
    real*8,allocatable :: syn_pr(:,:),ob1_pr(:,:)
    real*8,allocatable :: syn_pz(:,:),ob1_pz(:,:)
    real*8,allocatable :: syn_sr(:,:),ob1_sr(:,:)
    real*8,allocatable :: syn_st(:,:),ob1_st(:,:)
    real*8,allocatable :: syn_sz(:,:),ob1_sz(:,:)
    
    real*8 :: cc_pr,cc_pz,cc_sr,cc_sz,cc_st
    integer :: is_pr,is_pz,is_sr,is_st,is_sz
    integer :: counter,ist
    
    real*8 :: moment(6)
    real*8 :: norm_ob,norm_syn
    real*8 :: corrcoef,misfit,norm
    
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
    
    ! generate synthetic seismograms and align with observation
    counter = 0
    sca_M = 0.0
    do ist=1,nst
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
            else
                call align_seis(ob(ist),g(ist),nst,moment,'P','r',dt,pre_t(1),leng(1),0,    &
                            is_pz,ob1_pr(:,ist),syn_pr(:,ist))
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
            else
                call align_seis(ob(ist),g(ist),nst,moment,'S','r',dt,pre_t(2),leng(2),0,    &
                            is_sz,ob1_sr(:,ist),syn_sr(:,ist))
            end if
        end if
        
        !print *, is_pr,is_pz,is_sr,is_st,is_sz
        
        if (ob(ist).comp(1) /= 0) then
            norm_ob = norm_ob+sum(ob1_pr(:,ist)**2)
            norm_syn = norm_syn+sum(syn_pr(:,ist)**2)
        end if
        
        if (ob(ist).comp(2) /= 0) then
            norm_ob = norm_ob+sum(ob1_pz(:,ist)**2)
            norm_syn = norm_syn+sum(syn_pz(:,ist)**2)
        end if
        
        if (ob(ist).comp(3) /= 0) then
            norm_ob = norm_ob+sum(ob1_sr(:,ist)**2)
            norm_syn = norm_syn+sum(syn_sr(:,ist)**2)
        end if
        
        if (ob(ist).comp(4) /= 0) then
            norm_ob = norm_ob+sum(ob1_st(:,ist)**2)
            norm_syn = norm_syn+sum(syn_st(:,ist)**2)
        end if
        
        if (ob(ist).comp(5) /= 0) then
            norm_ob = norm_ob+sum(ob1_sz(:,ist)**2)
            norm_syn = norm_syn+sum(syn_sz(:,ist)**2)
        end if
        
        sca_M = sca_M + sqrt(norm_ob)/sqrt(norm_syn)
        
    end do    
    
    sca_M = sca_M/counter
    
    ! compute the misfit
    err = 0.0
    counter = 0
    do ist=1,nst
        !print *, any(ob(ist).comp /=0)
        if(.not. any(ob(ist).comp /=0)) cycle
        
        if (ob(ist).comp(1) /= 0) then ! Pr
            cc_pr = corrcoef(ob1_pr(:,ist),syn_pr(:,ist),leng(1))
            
            err = err + misfit(ob1_pr(:,ist),syn_pr(:,ist)*sca_M,leng(1),cc_pr,ob(ist).comp(1))
            
            counter = counter+1
        end if
        
        if (ob(ist).comp(2) /= 0) then ! Pz
            cc_pz = corrcoef(ob1_pz(:,ist),syn_pz(:,ist),leng(1))
            
            err = err + misfit(ob1_pz(:,ist),syn_pz(:,ist)*sca_M,leng(1),cc_pz,ob(ist).comp(2))
            
            counter = counter+1
        end if
    
        if (ob(ist).comp(3) /= 0) then ! Sr
            cc_sr = corrcoef(ob1_sr(:,ist),syn_sr(:,ist),leng(2))
            
            err = err + misfit(ob1_sr(:,ist),syn_sr(:,ist)*sca_M,leng(2),cc_sr,ob(ist).comp(3))
            !print *, cc_sr, misfit(ob1_sr(:,ist),syn_sr(:,ist)*sca_M,leng(2),cc_sr,ob(ist).comp(3)),err
            counter = counter+1
        end if
        
        if (ob(ist).comp(4) /= 0) then ! St
            cc_st = corrcoef(ob1_st(:,ist),syn_st(:,ist),leng(2))
            
            err = err + misfit(ob1_st(:,ist),syn_st(:,ist)*sca_M,leng(2),cc_st,ob(ist).comp(4))
            
            counter = counter+1
        end if
        
        if (ob(ist).comp(5) /= 0) then ! Sz
            cc_sz = corrcoef(ob1_sz(:,ist),syn_sz(:,ist),leng(2))
            
            err = err + misfit(ob1_sz(:,ist),syn_sz(:,ist)*sca_M,leng(2),cc_sz,ob(ist).comp(5))
            
            counter = counter+1
        end if
    end do
    err = err/counter
    
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
    end subroutine fitness
    
    function norm(x,nst)
    implicit none
    
    integer :: nst
    real*8 :: x(nst,1),tmp(1,1)
    
    real*8 :: norm
    
    tmp = matmul(transpose(x),x)
    norm = sqrt(tmp(1,1))
    
    return
    end function norm
    
    function misfit(ob,syn,leng,cc,weight)
    implicit none
    
    integer,intent(in) :: leng
    real*8,intent(in) :: ob(leng),syn(leng)
    real*8,intent(in) :: cc,weight
    
    real*8 :: misfit
    
    real*8 :: maxamp
    real*8 :: norm
    
    maxamp = max(maxval(abs(ob)),maxval(abs(syn)))
    misfit = norm(ob/maxamp-syn/maxamp,leng)*(1.0-cc)*weight
    !print *, misfit
    return
    end function misfit
    
    