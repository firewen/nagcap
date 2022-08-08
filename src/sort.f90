    subroutine sort(a,idx,n)
    implicit none
    integer :: n
    real*8 :: a(n)
    integer :: idx(n)
    
    integer :: i
    
    do i=1,n
        idx(i) = i
    end do
    
    call quick_sort(a,idx,n,1,n)
     
    return
    end subroutine sort
    
    recursive subroutine quick_sort(a,idx,n,s,e)
    implicit none
    integer :: n ! 表示类型的大小
    real*8 :: a(n) ! 存放数据的类型
    integer :: idx(n)
    integer :: s ! 传入的参数, 这一组的类型起始位置
    integer :: e ! 传入的参数, 这一组的类型结束位置
    
    integer :: l,r ! 用来找a(l)>k及a(r)<k时用的
    real*8 :: k ! 记录键值a(s)
    real*8 :: temp ! 交换两个数值时用的
    integer :: itmp
    
    ! 首先要先给定l,r的初值. l要从头开始,e则要从尾开始
    l=s
    r=e+1
    
    !print *,l,r
    ! right值 > left值 时才有必要进行排序
    if ( r<=l ) return

    k=a(s) ! 设定键值
    do while(.true.)
    ! 找出a(l)<k的所在
        do while( .true. )
            l=l+1
            if (l >= e) exit
            if (a(l) > k) exit
            !if ( (l>=e) .or. (a(l) > k) ) exit
        end do
    ! 找出a(r)>k的所在
        do while( .true. )
            r=r-1
            if (r <= s) exit
            if (a(r) < k) exit
            !if ( (r<=s) .or. (a(r) < k) ) exit
        end do
    ! 如果right 跑到 left的左边时, 循环就该结束了
        if ( r <= l ) exit
    ! 交换a(l),a(r)的数值
        temp=a(l)
        a(l)=a(r)
        a(r)=temp
        itmp = idx(l)
        idx(l) = idx(r)
        idx(r) = itmp
    end do
    ! 交换a(s),a(r)的数值
    temp=a(s)
    a(s)=a(r)
    a(r)=temp
    itmp = idx(s)
    idx(s) = idx(r)
    idx(r) = itmp
    ! 把r之前的数据重新分组,再做排序
    call quick_sort(a,idx,n,s,r-1)
    ! 把r之后的数据重新分组,再做排序
    call quick_sort(a,idx,n,r+1,e)
    return
    end subroutine quick_sort