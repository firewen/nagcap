    subroutine conv(x,h,y,n,m,l,ierror)
! --------------------------------------------------------------------
! Routine CONV: To implement Linear Convolution y(n)=x(n)*h(n)
! input parameters:
! x(n):L dimensioned real array,signal data is stored in x(0) to x(n-1).
! h(n):L dimensioned real array,impulse response is stored in h(0) to h(m-1).
! n   : the data length of x.
! m   : the data length of h.
! L   : the data length of y, L must be >=n+m-1
! output parameters:
! y(n):L dimensioned real array, y(n)=x(n)*h(n),n=0,...L-1.
!                                      in Chapter 1
!--------------------------------------------------------------------
    implicit none
    
    integer :: l,m,n,ierror
    real*8 :: x(0:n-1),h(0:m-1),xx(0:l-1),hh(0:l-1),y(0:l-1)
    
    integer :: i,k
    real*8 :: asum
    
    ierror = 0
    if (l < m+n-1) then
        ierror = 1
        return
    end if
    
    xx = 0.0
    xx(0:n-1) = x
    hh = 0.0
    hh(0:m-1) = h
    
    do k=0,l-1
        asum = 0.0
        do i=0,k
            asum = asum+xx(i)*hh(k-i)
        end do
        y(k) = asum
    end do
    
    return
    end subroutine conv