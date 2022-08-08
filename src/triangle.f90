    function triangle(lenstf,shift)
    implicit none
    
    integer :: lenstf,shift
    real*8 :: triangle(shift+lenstf+1)
    
    integer :: i
    real*8 :: a
    
    triangle = 0.0
    
    a = 2.0/lenstf
    do i=1,lenstf
        if (i <= lenstf/2) then
            triangle(shift+i) = a/(lenstf/2.0)*(i-1)
        else
            triangle(shift+i) = a/(lenstf/2.0)*(lenstf-i+1)
        end if
    end do
    
    return
    end function triangle