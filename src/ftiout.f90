    subroutine fitout(b,a,lb,la,x,n,y,ierror)
!-----------------------------------------------------------------------
! Routine FITOUT: To obtain the output of an IIR or a FIR filter;
!                 System function H(Z)=B(Z)/A(Z) has form as:
!                    b(0)+b(1)z^(-1)+ ... + b(lb)z^(-lb)
!              H(z)=------------------------------------
!                    1.0 +a(1)z^(-1)+ ... + a(la)z^(-la)
! input parameters:
! b :lb+1 dimensioned real array. b(0) to b(lb) is the coef. of B(z)
! a :la   dimensioned real array. a(1) to a(la) is the coef. of A(z)
! n :data point number of input and output signal.
! x :n dimensioned real array, x(0) to x(n-1):input signal ;
! output parameters:
! y :n dimensioned real array, y(0) to y(n-1):output signal;
!   if ierror=0: no errors detected,
!      ierror=1: output exceeds 1.E10.
! Note:
! If LA=0: For FIR system.
! If x(0)=1,x(n)=0 for n>0,the output is the impulse response h(n).
!                                          in Chapter 2
!-----------------------------------------------------------------------
    implicit none
    integer,intent(in) :: la,lb,n
    real*8,intent(in) :: x(0:n-1),b(0:lb),a(1:la)
    real*8,intent(out) :: y(0:n-1)
    integer,intent(out) :: ierror
    
    integer :: i,k,ki
    real*8 :: sum
    
    ierror = 1
    do k=0,n-1
        sum = 0.0
        do i=0,lb
            ki = k-i
            if (ki < 0) then
                exit
            end if
            if ((b(i) == 0.0) .or. (x(ki) == 0.0)) then
                cycle
            end if
            sum = sum+b(i)*x(ki)
            if(sum > 1.0e10) return
        end do
        if (la == 0) then
            y(k) = sum
        else
            do i=1,la
                ki = k-i
                if (ki < 0) then
                    exit
                end if
                if (a(i) == 0.0) then
                    cycle
                end if
                sum = sum-a(i)*y(ki)
                if(sum > 1.0e10) return
            end do
            y(k) = sum
        end if
    end do
    ierror = 0
    
    return
    end subroutine fitout

    