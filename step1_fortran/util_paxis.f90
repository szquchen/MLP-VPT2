module util_paxis

  contains
  !==================================================
  ! Diagonalize the matrix and return the eigenvalues and eigenvectors.
  ! The original matrix will be replaced with eigenvectors.
  !==================================================
  subroutine diag(dim,mat,w)
    integer,intent(in)::dim
    real,dimension(1:dim,1:dim),intent(inout)::mat
    real,dimension(1:dim),intent(out)::w
    ! ::::::::::::::::::::
    real,dimension(:),allocatable::work
    integer::lwork,info,i,j
  
    lwork=dim*dim*10
    allocate(work(1:lwork))
  
    call dsyev('v','u',dim,mat,dim,w,work,lwork,info)
  
    return
  end subroutine diag
  
  !==================================================
  ! Principal Axis Coordintates: X Y Z
  !==================================================
  subroutine prin_axis_coor(x,mass,eig)
    real,dimension(:),intent(inout)::x
    real,dimension(:),intent(in)::mass
    real,dimension(3),intent(out)::eig
    !::::::::::::::::::::
    real,dimension(3,size(mass))::xx
    real,dimension(size(mass),3)::xt
    real,dimension(3,3)::imat
    integer::i,natm
    real::tr

    natm = size(mass)
    do i=1,natm
       xx(:,i) = x(3*i-2:3*i)
    end do

    call com_coor(xx, mass)
  
    xt=transpose(xx)
    do i=1,natm
       xt(i,:)=xt(i,:)*mass(i)
    end do
    
    imat=-matmul(xx,xt)
    tr=imat(1,1)+imat(2,2)+imat(3,3)
    do i=1,3
       imat(i,i)=imat(i,i)-tr
    end do
  
    call diag(3,imat,eig)
  
    xx=matmul(transpose(imat),xx)
    do i=1,natm
       x(3*i-2:3*i) = xx(:,i)
    end do

    return
  end subroutine prin_axis_coor

  !==================================================
  ! Convert to Center of Mass Coordinates: X Y Z
  !==================================================
  subroutine com_coor(xx,mass)
    real,dimension(:,:),intent(inout)::xx
    real,dimension(:),intent(in)::mass
    !::::::::::::::::::::
    real,dimension(1:3)::xcm
    integer::i,natm
    real::m

    natm = size(mass)
    m=sum(mass)
    xcm=matmul(xx,mass)/m
  
    do i=1,natm
       xx(:,i)=xx(:,i)-xcm
    end do
  
    return
  end subroutine com_coor

end module util_paxis
