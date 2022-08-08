    subroutine findpeaks(src,lag,n,tshift,idx_out)
    integer :: n,tshift
    real*8 :: src(n)
    integer :: lag(n)
    
    integer :: sig(n-1)
    real*8 :: diff
    integer :: max_idx,idx_out
    integer,allocatable :: indmax(:),idx(:)
    real*8,allocatable :: srcmax(:)
    
    do i=1,n-1
        diff = src(i+1)-src(i)
        if (diff > 0) then
            sig(i) = 1
        elseif (diff < 0) then
            sig(i) = -1
        else
            sig(i) = 0
        end if
    end do
    
    max_idx = 0
    do i=1,n-2
        diff = sig(i+1)-sig(i)
        if (diff < 0) then
            max_idx = max_idx+1
        end if
    end do
    allocate(indmax(max_idx))
    allocate(idx(max_idx))
    allocate(srcmax(max_idx))
    max_idx = 0
    do i=1,n-2
        diff = sig(i+1)-sig(i)
        if (diff < 0) then
            max_idx = max_idx+1
            indmax(max_idx) = i+1
            srcmax(max_idx) = src(i+1)
        end if
    end do
    
    call sort(srcmax,idx,max_idx)
    
    idx_out = int(n/2)+1
    
    do i=max_idx,1,-1
        !print *, lag(indmax(idx(max_idx))),lag(indmax(idx(i)))
        if (lag(indmax(idx(max_idx))) < 0) then
            if ((lag(indmax(idx(i))) <0) .and. (lag(indmax(idx(i))) >-tshift)) then
                idx_out = indmax(idx(i))
                exit
            end if
        else 
            if ((lag(indmax(idx(i))) >=0) .and. (lag(indmax(idx(i))) <tshift)) then
                idx_out = indmax(idx(i))
                exit
            end if
        end if    
        !if ((lag(indmax(idx(i))) <tshift) .or. (lag(indmax(idx(i))) <tshift)) then
        !    idx_out = indmax(idx(i))
        !    exit
        !end if
    end do
    
    deallocate(indmax)
    deallocate(idx)
    deallocate(srcmax)
    return
    end subroutine findpeaks
    
        