    subroutine nmtensor(beta,gamma,str,dip,rake,tensor)
    ! beta : 0~pi
    ! gamma : -pi/6:pi/6
    
    implicit none
    
    real*8 :: beta, gamma, str, dip, rake
    real*8 :: tensor(6)
    
    real*8 :: n(3), v(3), b(3)
    real*8 :: diso(3,3), ddc(3,3), dclvd(3,3), eye(3,3), dd(3,3)
    real*8 :: iso, clvd
    
    integer :: i
    
    real*8 :: k_product(3,3)
    
    ! fault normal
    n(1) = -sin(dip)*sin(str)
    n(2) =  sin(dip)*cos(str)
    n(3) = -cos(dip)
    ! slip
    v(1) = cos(rake)*cos(str) + cos(dip)*sin(rake)*sin(str)
    v(2) = cos(rake)*sin(str) - cos(dip)*sin(rake)*cos(str)
    v(3) = -sin(rake)*sin(dip)
    ! the null vector
    b(1) = n(2)*v(3) - n(3)*v(2)
    b(2) = -(n(1)*v(3) - n(3)*v(1))
    b(3) = n(1)*v(2) - n(2)*v(1)
    
    eye = 0.0
    do i=1,3
        eye(i,i) = 1.0
    end do
    
    diso = sqrt(2.0/3.0)*eye
    ddc = matmul(reshape(n,(/3,1/)), reshape(v,(/1,3/))) + matmul(reshape(v,(/3,1/)), reshape(n,(/1,3/)))
    dclvd = 1.0/sqrt(3.0)*(2*matmul(reshape(b,(/3,1/)), reshape(b,(/1,3/))) &
            -matmul(reshape(v,(/3,1/)), reshape(v,(/1,3/)))-matmul(reshape(n,(/3,1/)), reshape(n,(/1,3/))))
    
    iso = cos(beta)
    clvd = sin(gamma)
    
    dd = iso*diso
    if ((1.0-iso*iso) > 0) then
        dd = dd + sqrt(1.0-iso*iso)*sqrt(1.0-clvd*clvd)*ddc
        if (abs(clvd) > 0.0001) then
            dd = dd + sqrt(1.0-iso*iso)*clvd*dclvd
        end if
    end if
    
    tensor(1) = dd(1,1)
    tensor(2) = dd(1,2)
    tensor(3) = dd(1,3)
    tensor(4) = dd(2,2)
    tensor(5) = dd(2,3)
    tensor(6) = dd(3,3)
    
    return
    end subroutine nmtensor
    
    
