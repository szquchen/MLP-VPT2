module pes_shell
  use constants
  use pes_basis
  implicit none

  real(kind=8)::coef(49978)


contains
  !==========================!
  ! read in the coefficients !
  !==========================!
  subroutine pes_init()
    character(len=90)::path
    integer::i

    path = "coeff_Aspirin.dat"
    open(20,file=trim(path),status="old")
    do i=1,size(coef)
       read(20,*) coef(i)
    end do
    close(20)

    return
  end subroutine

  !================================!
  ! calculate the potential energy !
  !   x(3,15): coordinates in bohr !
  !================================!
  function pot(x)
    real(kind=8)::x(63),xyz(21,3)
    real(kind=8)::m(11247), p(49978), morse(210)
    real(kind=8)::pot
    integer::i

    do i=1,21!end of do should be # of atoms
       xyz(i,:) = x(3*i-2:3*i)
    end do

    call get_x(xyz, morse)
    call evmono(morse, m)
    call evpoly(m, p)
    pot = dot_product(coef, p)

    return
  end function

  !================================!
  ! calculate the gradients        !
  !   x(3,15): coordinates in bohr !
  !================================!
  function gradient(x)
    real(kind=8)::x(63), xyz(21,3), gradient(63)
    real(kind=8)::m(11247), p(49978), morse(210)
    integer::i

    do i=1,21
       xyz(i,:) = x(3*i-2:3*i)
    end do

    call get_x(xyz, morse)
    call evmono(morse, m)
    call evpoly(m, p)
    call derivative_reverse(coef,m,p,xyz,gradient)

    return
  end function

  !=============================!
  ! Hessian using the gradient  !
  !=============================!
  subroutine hessian(x,H)
    real,dimension(:),intent(in)::x
    real,dimension(:,:),intent(inout)::H
    !:::::::::::::::::::::::::::
    real,dimension(1:size(x))::xt,g_1f,g_1b,g_2f,g_2b
    real,parameter::eps=0.001
    integer::ndim,i

    ndim = size(x)
    do i = 1, ndim
       xt=x;  xt(i)=xt(i)+2.0*eps;  g_2f=gradient(xt)
       xt=x;  xt(i)=xt(i)+1.0*eps;  g_1f=gradient(xt)
       xt=x;  xt(i)=xt(i)-1.0*eps;  g_1b=gradient(xt)
       xt=x;  xt(i)=xt(i)-2.0*eps;  g_2b=gradient(xt)

       H(:,i) = (-g_2f + 8.0*g_1f - 8.0*g_1b + g_2b) / (12.0*eps)
    end do

    return
  end subroutine hessian

end module
