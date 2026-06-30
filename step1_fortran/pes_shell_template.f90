module pes_shell
  use constants
  implicit none

contains
  !=========================================!
  ! Use this subroutine to initialize the   !
  ! PES (e.g. read in parameters, etc.)     !
  !=========================================!
  subroutine pes_init()

    return
  end subroutine

  !=========================================!
  ! Calculate the potential energy given    !
  ! the Cartesian coordinates               !
  !   input x(1:3*Natom): Cartesian in bohr !
  !   output pot: energy in hartree         !
  !=========================================!
  function pot(x)

    ! Your code here

    return
  end function

  !=========================================!
  ! Calculate gradient (Not force) given    !
  ! Cartesian coordinates                   !
  !   input: x(1:3*Natom) Cartesian in bohr !
  !   output: gradient(1:3*Natom) gradient  !
  !           in atomic unit (hartree/bohr) !
  !=========================================!
  function gradient(x)

    return
  end function

  !=========================================!
  ! Calculating Hessian matrix using finite !
  ! difference of gradient                  !
  !   input: x(1:3*Natom) Cartesian in bohr !
  !   output: H(3*Natom, 3*Natom) Hessian   !
  !           in hartree/(bohr^2)           !
  !=========================================!
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
