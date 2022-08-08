    function detrend(y,n)
    implicit none
    
    integer :: n
    real*8 :: y(n)
    real*8 :: x(n)
        
    real*8 :: s,t
    real*8 :: a,b,c,d,e,f
    integer :: i
    
    real*8 :: detrend(n)
    
    do i=1,n
        x(i) = real(i,kind=8)
    end do
    
    a = 0; b = 0; e = 0; f = 0
    do i=1,n
        a = a+x(i)*x(i)
        b = b+x(i)
        e = e+x(i)*y(i)
        f = f+y(i)
    end do
    c = b
    d = n
    
    s = (b*f-e*d)/(b*c-a*d)
    t = (e*c-a*f)/(b*c-a*d)
    detrend = y-(s*x+t)
    
    return
    end function detrend