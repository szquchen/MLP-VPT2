!=============================================================================!
! this code is used to do geometry optimization using Newton's Method.        !
! Zhen Xie@Bowman Group                                                       !
! zxie3@emory.edu                                                             !
!                                                                             !
! Dec. 16 2004                                                                !
! Emerson Center , Emory University                                           !
!=============================================================================!
module util_opt
  use pes_shell
  implicit none
   
  contains
  !================================================================!
  ! print the molecule to a file using a coordinates vector        !
  ! the output file includes the geometry and enrgy in xyz         !
  ! format                                                         !
  !================================================================!
  subroutine prtxmol(x,symbs,file)
    double precision,dimension(:),intent(in)::x
    character(len=2),dimension(:),intent(in)::symbs
    integer::file
    ! ::::::::::::::::::::
    integer::i,natm
    
    natm=size(x,1)/3
    
    write(file,'(I2)') natm
    write(file,*) "Angstroms"
    
    do i=1,natm
       write(file,'(A,3F20.14)') symbs(i),x(3*i-2:3*i)*auang
    end do

    return
  end subroutine prtxmol

  !===============================================!
  ! optimize the geometry using the newton        !
  ! method p=p-LD^{-1}L'g                         !
  !===============================================!
  subroutine optg(p)
    double precision,dimension(:),intent(inout)::p
    ! ::::::::::::::::::::
    double precision,dimension(1:size(p))::grad,w,disp,tp,tdisp
    double precision,dimension(1:size(p),1:size(p))::H,invw
    double precision::maxd,enrg,maxg,imaxg,lambda
    integer::n,i

    n=0
    grad = gradient(p)       !calculate the gradient vector
    maxg=maxval(abs(grad))
    imaxg=maxg

    enrg=pot(p)

    write(*,'(A)') "1. Start geometry optimization..."
    write(*,*)
    write(*,'(A,I3,A,F10.7,A,F15.6)') &
         " STEP:",n," MAX_GRAD:",maxg," ENERGY:", enrg
    
    do while(maxg>1.0e-6 .and. n<=100)
       call hessian(p,H)        !calculate the hessian Matrix
       call diag_hessian(H,w)   !diagnolize the hessian Matrix
       call inverse_w(w,invw) !get the inverse matrix of eigenvalues
       call displacement(grad,H,invw,disp) !calculae the displacement 
       
       lambda=1.0               !step length in p=p-\lambda*LD^{-1}L'g
       do while(lambda>1.0e-5)
          tdisp=lambda*disp
          tp=p                  
          call updatep(tp,tdisp) !update the point
          grad = gradient(tp) !calculate the gradient vector
          maxg=maxval(abs(grad)) !maximum gradient
          if(maxg>imaxg) then
             lambda=lambda*0.5
             cycle
          else
             disp=tdisp
             p=tp
             imaxg=maxg
             exit
          end if
       end do
       
       if(lambda<=1.0e-4) then
          write(*,'(A)') "Warning: Optimization not converged!"
          return
       end if
       maxd=maxval(abs(disp)) !maximum displacement
       enrg=pot(p)
       n=n+1

       write(*,'(A,I3,A,F10.7,A,F10.7,A,F12.6)') &
            "STEP:",n," MAX_GRAD:",imaxg," MAX_DISP:", maxd," ENERGY:", enrg
       
    end do

    write(*,*)
    write(*,'(A)') "Optimization Completed."
    
    return
  end subroutine optg

  !=====================================================!
  ! calculate the inverse of the eigen values w         !
  ! if the w(i) is too small, set it to 0.              !
  ! in the real case, there should be 6 zeros in        !
  ! the w(i)s. The return should be a diagnonal         !
  ! matrix                                              !
  !=====================================================!
  subroutine inverse_w(w,invw)
    double precision,dimension(:),intent(in)::w
    double precision,dimension(:,:),intent(out)::invw
    ! ::::::::::::::::::::
    integer::i,dim
   
    dim=size(w,1)
    invw=0
    
    do i=7,dim
       invw(i,i)=1/w(i)
    end do
    
    return
  end subroutine inverse_w

  !=====================================================!
  ! update the current point p to the next point        !
  ! p=p+disp                                            !
  !=====================================================!
  subroutine updatep(p,disp)
    double precision,dimension(:),intent(inout)::p
    double precision,dimension(:),intent(in)::disp
    ! ::::::::::::::::::::
    p=p+disp
    return
  end subroutine updatep

  !==========================================================!
  ! calculate the displacement of point p: disp using        !
  ! the diagnolized hessian                                  !
  ! disp=-L*invw*L'g                                         !  
  !==========================================================!
  subroutine displacement(grad,norms,invw,disp)
    double precision,dimension(:),intent(in)::grad
    double precision,dimension(:,:),intent(in)::norms,invw
    double precision,dimension(:),intent(out)::disp
    ! ::::::::::::::::::::
    
    disp=matmul(transpose(norms),grad)
    disp=matmul(invw,disp)
    disp=matmul(norms,disp)

    disp=-disp

    return
  end subroutine displacement

  !==================================================!
  ! diagonalize the hessian matrix and return        !
  ! the eigen value in freq and eigenvectors         !
  ! in norms. The original Hesion matrix             !
  ! will be destroied                                !
  !==================================================!
  subroutine diag_hessian(H,w)
    double precision,dimension(:,:),intent(inout)::H
    double precision,dimension(:),intent(out)::w
    ! ::::::::::::::::::::
    double precision,dimension(:),allocatable::work
    integer::dim,lwork,info,i,j
    
    dim=size(H,1)
    lwork=dim*dim*10;
    allocate(work(1:lwork))
    
    call dsyev('v','u',dim,H,dim,w,work,lwork,info) 
    
    return
  end subroutine diag_hessian

end module util_opt
