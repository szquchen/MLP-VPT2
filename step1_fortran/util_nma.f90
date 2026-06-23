module util_nma
use pes_shell

  real,dimension(:),allocatable::x0,mass,freq,q
  real,dimension(:,:),allocatable::vec,Hc
  character(len=2),dimension(:),allocatable::symb

contains
  !==========================================
  ! calculate the mass weighted hessian at p
  !==========================================
  subroutine mw_hessian(p,H)
    real,dimension(:),intent(in)::p
    real,dimension(:,:),intent(inout)::H
    !:::::::::::::::::::
    real::rmass
    integer::dim,i,j

    dim=size(p)
    call hessian(p,H)

    do i=1,dim
      do j=1,dim
        rmass=sqrt(mass(ceiling(i/3.0)))*sqrt(mass(ceiling(j/3.0)))
        H(i,j)=H(i,j)/rmass
      end do
    end do

    return
  end subroutine mw_hessian

  !==================================================
  ! diagonalize the hessian matrix and return
  ! the eigen value and eigenvectors.
  ! The original Hessian matrix  will be destroied
  !==================================================
    subroutine diag_mwhessian(H,w)
    real,dimension(:,:),intent(inout)::H
    real,dimension(:),intent(out)::w
    ! ::::::::::::::::::::
    real,dimension(:),allocatable::work
    integer::dim,lwork,info,i,j

    dim=size(H,1)
    lwork=dim*dim*10
    allocate(work(1:lwork))

    call dsyev('v','u',dim,H,dim,w,work,lwork,info)

    do i=1,dim
       w(i)=sign(sqrt(abs(w(i)))*aucm,w(i))
    end do

    return
  end subroutine diag_mwhessian

  subroutine q2c(qq,xx)
    !==============================================!
    ! convert the normal coordinates to the        !
    ! usual cardesian coordinates                  !
    !==============================================!
    real,dimension(:),intent(in)::qq !Normal Coordinates
    real,dimension(:),intent(out)::xx ! Cartesian
    ! ::::::::::::::::::::
    real,dimension(1:size(xx))::mw !mass weighted displacements

    call q2mw(qq,mw)
    call mw2c(mw,xx)

    return
  end subroutine q2c

  subroutine q2mw(qq,mw)
    !=============================================!
    ! convert the normal coordinate to the        !
    ! mass weighted coordinate                    !
    !=============================================!
    real,dimension(:),intent(in)::qq ! normal coordinate
    real,dimension(:),intent(out)::mw
    ! ::::::::::::::::::::

    mw=matmul(vec,qq)

    return
  end subroutine q2mw

  subroutine mw2c(mw,xx)
    !=================================================!
    ! convert the mass weighted coordinates to        !
    ! the cardesian coordinates                       !
    !=================================================!
    real,dimension(:),intent(in)::mw
    real,dimension(:),intent(out)::xx
    ! ::::::::::::::::::::
    integer::i

    do i=1,size(xx)
       xx(i)=mw(i)/sqrt(mass(ceiling(i/3.0)))
    end do
    xx=xx+x0

    return
  end subroutine mw2c

  !==========================================!
  ! Given the normal coordinates, calculate  !
  ! the potential energy                     !
  !==========================================!
  function fq(qq)
    real,dimension(:),intent(in)::qq
    !:::::::::::::::::::::::::::::
    real,dimension(size(x0))::xx
    integer::natm,i
    real::fq

    call q2c(qq,xx) ! convert normal coord to Cartesian

    fq = pot(xx) ! replace with your potential routine

    return
  end function

  !=============================================!
  ! transform the cartesian grad to normal grad !
  !=============================================!
  function gq(qq)
    real,dimension(:),intent(in)::qq ! normal coordinate
    !::::::::::::::::::::::::::::
    real,dimension(size(x0))::xx,gx
    integer::natm,i
    real,dimension(size(qq))::gq

    call q2c(qq,xx) ! convert normal coord to Cartesian

    gx = gradient(xx) ! calculate grad in Cartesian space

    ! Then transform gx to grad in normal coordinate space
    do i=1,size(gx)
       gx(i) = gx(i) / sqrt(mass(ceiling(i/3.0)))
    end do
    gq = matmul(transpose(vec), gx)

    return
  end function

  !============================================!
  ! Compute Coriolis coupling constants zeta   !
  ! given mass-weighted normal mode vectors    !
  ! vec(3*natom, nmode)                        !
  !============================================!
  subroutine compute_coriolis(vec, zeta)
    real,dimension(:,:),intent(in)::vec
    real,dimension(3,size(vec,2),size(vec,2))::zeta
    !::::::::::::::::::
    integer::iat, ix, iy, iz, k, l, natm, nmode

    natm = size(vec,1) / 3
    nmode = size(vec,2)

    do iat = 1, natm

      ix = 3*(iat-1) + 1
      iy = ix + 1
      iz = ix + 2

      do k=1,nmode
        do l=1,nmode

          ! zeta_x
          zeta(1,k,l) = zeta(1,k,l) + vec(iy,k)*vec(iz,l) - vec(iz,k)*vec(iy,l)

          ! zeta_y
          zeta(2,k,l) = zeta(2,k,l) + vec(iz,k)*vec(ix,l) - vec(ix,k)*vec(iz,l)

          ! zeta_z
          zeta(3,k,l) = zeta(3,k,l) + vec(ix,k)*vec(iy,l) - vec(iy,k)*vec(ix,l)

        end do
      end do
    end do

    return
  end subroutine

end module util_nma
