    subroutine get_par(filename,ddir,gdir,leng,flh_p,flh_s,weight,lenstf,evdp,tshift,lb,hb)
    ! Input: 
    ! filename  : the file name wich includes all inversion parameters
    ! Output:
    ! ddir      : observation data directory
    ! gdir      : Green's functions directory
    ! leng      : the length for P and S phases
    ! flh_p     : the cutoff corner frequency of a filter for P
    ! flh_s     : the cutoff corner frequency of a filter for S
    ! weight    : a file path which includes station list and weight for each station
    ! lenstf    : the length of stf
    ! evdp      : the event depth (start end dh)
    ! tshift    : the permitted time shift for P and S phases
    ! lb        : the lower boundary for strike dip rake beta and gamma
    ! hb        : the upper boundary for strike dip rake beta and gamma
    implicit none
    
    character(len=*) :: filename,ddir,gdir,weight
    integer :: leng(2),lenstf,tshift(2)
    real*8 :: flh_p(2),flh_s(2),evdp(3),lb(5),hb(5)
    
    character*80 :: list,tmp
    
    open(10,file=filename,form='formatted',status='old')
    read(10,'(a80)') list
    read(10,'("ddir = "a30)') ddir
    read(10,'(a80)') list
    read(10,'("gdir = "a30)') gdir
    read(10,'(a80)') list
    read(10,'("len = "a80)') tmp
    read(tmp,*) leng
    read(10,'(a80)') list
    read(10,'("flh_p = "a80)') tmp
    read(tmp,*) flh_p
    read(10,'(a80)') list
    read(10,'("flh_s = "a80)') tmp
    read(tmp,*) flh_s
    read(10,'(a80)') list
    read(10,'("weight = "a40)') weight
    read(10,'(a80)') list
    read(10,'("lenstf = "I)') lenstf
    read(10,'(a80)') list
    read(10,'("evdp = "a80)') tmp
    read(tmp,*) evdp
    read(10,'(a80)') list
    read(10,'("tshift = "a80)') tmp
    read(tmp,*) tshift
    read(10,'(a80)') list
    read(10,'("lb = "a80)') tmp
    read(tmp,*) lb
    read(10,'(a80)') list
    read(10,'("hb = "a80)') tmp
    read(tmp,*) hb
    
    close(10)
    
    return
    end subroutine get_par