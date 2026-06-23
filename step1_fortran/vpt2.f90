program main
use util_opt
use util_paxis
use util_force

  implicit none
  real,dimension(:),allocatable :: fc1
  real,dimension(:,:),allocatable :: fc2, chi, H0
  real,dimension(:,:,:),allocatable :: fc3, Hp1, Hp2, Hm1, Hm2
  real,dimension(:,:),allocatable :: fc4 ! only need Uiiii and Uiijj
  real,dimension(:,:,:),allocatable :: zeta ! Coriolis coupling constants
  real :: dijk, wijk, I_pa(3), rotc(3), zpe, hzpe, chi0, va, temp
  integer :: natm, n, i, j, k, l
  character(len=32) :: fname, molname

  call getarg(1,fname) ! Command line argument: input file name
  i = index(fname, '.', .true.)
  molname = fname(1:i-1)

  write(*,'(A)') "***************************************************"
  write(*,'(A)') "Start VPT2 calculation"
  write(*,'(A)') "***************************************************"
  write(*,*)
  call pes_init()

  open(10,file=trim(fname),status='old')
  ! output files:
  ! molname_opt.out: optimized geometry in principal axis frame
  ! molname_nma.out: normal mode frequencies and eigenvectors
  ! molname_cubics.out: cubic force constants
  ! molname_quartics.out: quartic force constants
  open(11,file=trim(molname)//'_opt.out',status='unknown')
  open(12,file=trim(molname)//'_nma.out',status='unknown')
  open(13,file=trim(molname)//'_cubics.out',status='unknown')
  open(14,file=trim(molname)//'_quartics.out',status='unknown')

  read(10,*) natm
  read(10,*)
  n = 3*natm - 6 ! number of normal modes
  
  allocate(x0(3*natm), mass(natm), symb(natm))
  allocate(Hc(3*natm,3*natm), freq(3*natm), vec(3*natm,n))
  allocate(q(n), fc1(n), fc2(n,n), fc3(n,n,n), fc4(n,n), zeta(3,n,n))
  allocate(chi(n,n))

  ! Hessians at displaced configs precomputed to avoid repeating these
  allocate(Hp1(n,n,n), Hp2(n,n,n), Hm1(n,n,n), Hm2(n,n,n), H0(n,n))

  ! read in the Cartesian coordinates
  do i=1,natm
     read(10,*) symb(i),x0(3*i-2:3*i)

     select case(symb(i))
     case("C")
        mass(i)=c_mass
     case("O")
        mass(i)=o_mass
     case("H")
        mass(i)=h_mass
     case("N")
        mass(i)=n_mass
     end select

     mass(i)=mass(i)*emass
  end do
  x0 = x0 / auang
  close(10)

  ! 1) Geometry optimization and rotate to principal axis frame
  call optg(x0)
  call prin_axis_coor(x0, mass, I_pa)
  call prtxmol(x0, symb, 11)
  write(*,'(A)') "Optimized geometry written to "//trim(molname)//"_opt.out"
  write(*,'(A)') "***************************************************"
  write(*,*)
  close(11)

  ! 2) Normal mode analysis
  call mw_hessian(x0, Hc)
  call diag_mwhessian(Hc, freq)
  write(*,'(A)') "2. Normal mode analysis..."
  write(*,*)

  do i=1,n
     vec(:,i) = Hc(:,i+6)
  end do
  
  ! Hessian in Q at minimum H0
  q = 0.0
  call second_order(q,H0)
  write(12,'(A)') "Frequencies (cm^-1):"
  do i=1,n
     fc2(i,i) = sqrt(H0(i,i))
     write(12,'(I5,F18.8)') i,fc2(i,i) * aucm
  end do
  write(12,*)
  do i=1,n
     write(12,'(A,I5)') "Mode", i
     do j=1,natm
        write(12,'(3F20.10)') vec(3*j-2:3*j,i)
     end do
  end do
  write(*,'(A)') "Done. Results written to "//trim(molname)//"_nma.out"
  write(*,'(A)') "***************************************************"
  write(*,*)
  close(12)

  ! Precompute the Hessians and save them
  do k=1,n
     q = 0.0
     q(k) = det34
     call second_order(q,Hp1(:,:,k))
     q(k) = 2.0*det34
     call second_order(q,Hp2(:,:,k))
     q(k) = -det34
     call second_order(q,Hm1(:,:,k))
     q(k) = -2.0*det34
     call second_order(q,Hm2(:,:,k))
  end do

  ! third-order force constants
  write(*,'(A)') "3. Computing cubic force constants"
  write(*,*)
  call third_order(Hp1,Hm1,Hp2,Hm2,fc3)
  write(13,'(A)') "#  i  j  k        Phi_ijk   (Hartree / Q^3)"
  do i=1,n
    do j=i,n
      do k=j,n
        ! symmetrize the constants for all 6 permutations of (i,j,k)
        temp = ( fc3(i,j,k) + fc3(i,k,j) + fc3(j,i,k) + &
                 fc3(j,k,i) + fc3(k,i,j) + fc3(k,j,i) ) / 6.0
        write(13,'(3I5,E24.12)') i,j,k,temp

        ! convert from Q to q, and Eh to cm-1
        temp = temp / sqrt(fc2(i,i)*fc2(j,j)*fc2(k,k)) * aucm

        fc3(i,j,k) = temp
        fc3(i,k,j) = temp
        fc3(j,i,k) = temp
        fc3(j,k,i) = temp
        fc3(k,i,j) = temp
        fc3(k,j,i) = temp
      end do
    end do
  end do
  write(*,'(A)') "Done. Results written to "//trim(molname)//"_cubics.out"
  write(*,'(A)') "***************************************************"
  write(*,*)
  close(13)

 ! fourth-order force constants
  write(*,'(A)') "4. Computing quartic force constants"
  call fourth_order(H0,Hp1,Hm1,Hp2,Hm2,fc4)
  write(14,'(A)') "#  i  j  k  l      Phi_ijkl   (Hartree / Q^4)"
  do i=1,n
    do j=i,n
      write(14,'(4I5,E24.12)') i,i,j,j,fc4(i,j)
      fc4(i,j) = fc4(i,j) / sqrt(fc2(i,i)**2*fc2(j,j)**2) * aucm
      fc4(j,i) = fc4(i,j)
    end do
  end do
  write(*,'(A)') "Done. Results written to "//trim(molname)//"_quartics.out"
  write(*,'(A)') "***************************************************"
  write(*,*)
  close(14)

  ! convert second order constants to cm-1
  do i=1,n
     fc2(i,i) = fc2(i,i) * aucm
  end do

  ! calculate the VPT2 energies without resonance treatment
  write(*,'(A)') "5. VPT2 fundamentals without resonance:"
  write(*,*)
  
  ! rotational constants and Coriolis coupling zeta
  call compute_coriolis(vec, zeta)
  do i=1,3
     rotc(i) = 0.5 / I_pa(i) * aucm
  end do
  write(*,'(A)') "Rotational constants (cm-1):"
  write(*,'(3F15.8)') rotc(1:3)
  write(*,*)

  chi0 = 0.0
  do i=1,n
    chi(i,i) = fc4(i,i)
    chi0 = chi0 + fc4(i,i) - 7.0 / 9.0 * fc3(i,i,i)*fc3(i,i,i) / fc2(i,i)
    do j=1,n
      chi(i,i) = chi(i,i) - (8.0*fc2(i,i)*fc2(i,i) - 3.0*fc2(j,j)*fc2(j,j)) * &
                 fc3(i,i,j)*fc3(i,i,j) / fc2(j,j) / &
                 (4.0*fc2(i,i)*fc2(i,i)-fc2(j,j)*fc2(j,j))
   
      if (j .ne. i) then
        chi(i,j) = fc4(i,j)
        chi0 = chi0 + 3.0 * fc2(i,i)*fc3(i,j,j)*fc3(i,j,j) / &
               (4.0*fc2(j,j)*fc2(j,j) - fc2(i,i)*fc2(i,i))

        do k=1,n
          dijk = (fc2(i,i)+fc2(j,j)-fc2(k,k))*(fc2(i,i)+fc2(j,j)+fc2(k,k))* &
               & (fc2(i,i)-fc2(j,j)+fc2(k,k))*(fc2(i,i)-fc2(j,j)-fc2(k,k))
          wijk = fc2(i,i)*fc2(i,i) + fc2(j,j)*fc2(j,j) - fc2(k,k)*fc2(k,k)
          chi(i,j) = chi(i,j) - fc3(i,i,k)*fc3(j,j,k)/fc2(k,k) + &
                   & 2.0*fc2(k,k)*wijk*fc3(i,j,k)*fc3(i,j,k) / dijk
        end do

        do k=1,3 ! Coriolis contribution
          chi(i,j) = chi(i,j) + 4.0*(fc2(i,i)**2+fc2(j,j)**2) / &
                     (fc2(i,i)*fc2(j,j))*rotc(k)*zeta(k,i,j)**2
        end do

      chi(i,j) = chi(i,j) / 4.0
      end if
    end do
    chi(i,i) = chi(i,i) / 16.0

    do j=i+1,n
      do k=j+1,n
        dijk = (fc2(i,i)+fc2(j,j)-fc2(k,k))*(fc2(i,i)+fc2(j,j)+fc2(k,k))* &
             & (fc2(i,i)-fc2(j,j)+fc2(k,k))*(fc2(i,i)-fc2(j,j)-fc2(k,k))
        chi0 = chi0 - 16.0 * fc2(i,i)*fc2(j,j)*fc2(k,k) * &
               fc3(i,j,k)*fc3(i,j,k) / dijk
      end do
    end do

  end do

  do k=1,3
    temp = 1.0
    do i=1,n-1
      do j=i+1,n
         temp = temp + 2.0*zeta(k,i,j)**2
      end do
    end do
    chi0 = chi0 - 16.0*rotc(k)*temp
  end do
  chi0 = chi0 / 64.0

  ! Zero-point energies
  zpe = chi0 ! anharmonic ZPE
  hzpe = 0.0 ! harmonic ZPE
  do i=1,n
     hzpe = hzpe + 0.5*fc2(i,i)

     zpe = zpe + 0.5*fc2(i,i)
     do j=i,n
        zpe = zpe + chi(i,j)*0.25
     end do
  end do
  write(*,'(A,F12.2)') "Harmonic ZPE (cm-1): ", hzpe
  write(*,'(A,F12.2)') "Anharmonic ZPE (cm-1): ", zpe
  write(*,'(A,F12.2)') "Shift (cm-1): ", zpe-hzpe
  write(*,*)

  ! The fundamental energies
  write(*,'(A)') "Fundamentals (cm-1):"
  write(*,'(A)') "# Mode  Harmonic    VPT2     Shift"
  write(*,'(A)') "------------------------------------"
  do i=1,n
     va = fc2(i,i) + 2.0*chi(i,i)
     do j=1,n
        if (j .ne. i) then
           va = va + 0.5 * chi(i,j)
        end if
     end do
     write(*,'(I5,3F10.2)') i,fc2(i,i),va,va-fc2(i,i)
  end do

  write(*,*)
  write(*,'(A)') "VPT2 done."
  write(*,'(A)') "For DVPT2 and GVPT2 treatment of resonance, copy"
  write(*,'(A)') "all *.out files to step2 and run the Python script."
  write(*,'(A)') "***************************************************"

  deallocate(x0, mass, symb)
  deallocate(Hc, freq, vec)
  deallocate(q, fc1, fc2, fc3, fc4, zeta)
  deallocate(Hp1, Hp2, Hm1, Hm2, H0, chi)
  
end program
