module util_force
use util_nma

  real,parameter::det2=0.20
  real,parameter::det34=0.30

contains
  !==============================!
  ! calculating first derivative !
  !==============================!
  subroutine first_order(q,f1)
    real,dimension(:),intent(in)::q
    real,dimension(size(q)),intent(out)::f1
    !::::::::::::::

    f1 = gq(q)

    return
  end subroutine

  !=================================!
  ! calculating second derivatives  !
  !=================================!
  subroutine second_order(q,f2)
    real,dimension(:),intent(in)::q
    real,dimension(size(q),size(q)),intent(out)::f2
    !::::::::::::::
    real,dimension(size(q))::qt,f2f,f1f,f1b,f2b
    integer::i,j

    do j=1,size(q)
       qt = q
       qt(j) = qt(j) + 2.0*det2
       call first_order(qt,f2f)

       qt = q
       qt(j) = qt(j) + det2
       call first_order(qt,f1f)

       qt = q
       qt(j) = qt(j) - 2.0*det2
       call first_order(qt,f2b)

       qt = q
       qt(j) = qt(j) - det2
       call first_order(qt,f1b)

       f2(:,j) = (-f2f + 8.0*f1f - 8.0*f1b + f2b) / (12.0*det2)
    end do

    do i=1,size(q)-1
       do j=i+1,size(q)
          temp = (f2(i,j) + f2(j,i)) / 2.0
          f2(i,j) = temp
          f2(j,i) = temp
       end do
    end do

    return
  end subroutine

  !===============================!
  ! calculating third derivatives !
  !===============================!
  subroutine third_order(Hp1,Hm1,Hp2,Hm2,f3)
    real,dimension(:,:,:),intent(out)::f3
    real,dimension(:,:,:),intent(in)::Hp1,Hm1,Hp2,Hm2
    !::::::::::::::
    real,dimension(size(q),size(q))::f2f,f1f,f1b,f2b
    integer::i,j,k
    real::temp

    do k=1,size(f3,3)
       f3(:,:,k) = (-Hp2(:,:,k) + 8.0*Hp1(:,:,k) &
                    -8.0*Hm1(:,:,k) + Hm2(:,:,k)) / (12.0*det34)
    end do

    return
  end subroutine

  !================================!
  ! calculating fourth derivatives !
  !================================!
  subroutine fourth_order(H0,Hp1,Hm1,Hp2,Hm2,f4)
    real,dimension(:,:),intent(out)::f4
    real,dimension(:,:),intent(in)::H0
    real,dimension(:,:,:),intent(in)::Hp1,Hm1,Hp2,Hm2
    !::::::::::::::
    integer::i,j
    real::val

    do i=1,size(f4,1)
       f4(i,i) = (-Hp2(i,i,i) + 16.0*Hp1(i,i,i) - 30.0*H0(i,i) &
                  +16.0*Hm1(i,i,i) - Hm2(i,i,i)) / (12.0*det34**2)
    end do

    do i=1,size(f4,1)-1
       do j=i+1,size(f4,1)
          val = (-Hp2(i,i,j) + 16.0*Hp1(i,i,j) - 30.0*H0(i,i) &
                 +16.0*Hm1(i,i,j) - Hm2(i,i,j)) / (12.0*det34**2)
          f4(i,j) = val
          f4(j,i) = val
       end do
    end do

    return
  end subroutine

end module util_force
