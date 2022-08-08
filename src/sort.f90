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
    integer :: n ! ��ʾ���͵Ĵ�С
    real*8 :: a(n) ! ������ݵ�����
    integer :: idx(n)
    integer :: s ! ����Ĳ���, ��һ���������ʼλ��
    integer :: e ! ����Ĳ���, ��һ������ͽ���λ��
    
    integer :: l,r ! ������a(l)>k��a(r)<kʱ�õ�
    real*8 :: k ! ��¼��ֵa(s)
    real*8 :: temp ! ����������ֵʱ�õ�
    integer :: itmp
    
    ! ����Ҫ�ȸ���l,r�ĳ�ֵ. lҪ��ͷ��ʼ,e��Ҫ��β��ʼ
    l=s
    r=e+1
    
    !print *,l,r
    ! rightֵ > leftֵ ʱ���б�Ҫ��������
    if ( r<=l ) return

    k=a(s) ! �趨��ֵ
    do while(.true.)
    ! �ҳ�a(l)<k������
        do while( .true. )
            l=l+1
            if (l >= e) exit
            if (a(l) > k) exit
            !if ( (l>=e) .or. (a(l) > k) ) exit
        end do
    ! �ҳ�a(r)>k������
        do while( .true. )
            r=r-1
            if (r <= s) exit
            if (a(r) < k) exit
            !if ( (r<=s) .or. (a(r) < k) ) exit
        end do
    ! ���right �ܵ� left�����ʱ, ѭ���͸ý�����
        if ( r <= l ) exit
    ! ����a(l),a(r)����ֵ
        temp=a(l)
        a(l)=a(r)
        a(r)=temp
        itmp = idx(l)
        idx(l) = idx(r)
        idx(r) = itmp
    end do
    ! ����a(s),a(r)����ֵ
    temp=a(s)
    a(s)=a(r)
    a(r)=temp
    itmp = idx(s)
    idx(s) = idx(r)
    idx(r) = itmp
    ! ��r֮ǰ���������·���,��������
    call quick_sort(a,idx,n,s,r-1)
    ! ��r֮����������·���,��������
    call quick_sort(a,idx,n,r+1,e)
    return
    end subroutine quick_sort