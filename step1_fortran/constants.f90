module constants

  ! The mass of the most abundant isotope of each element
  real(kind=8),parameter :: c_mass= 12.0000000 ! 12C
  real(kind=8),parameter :: h_mass=  1.0078250 ! H
  real(kind=8),parameter :: d_mass=  2.0141018 ! D
  real(kind=8),parameter :: o_mass= 15.9949146 ! 16O
  real(kind=8),parameter :: n_mass= 14.003074  ! 14N
  real(kind=8),parameter :: f_mass= 18.9984032 ! 19F
  real(kind=8),parameter ::si_mass= 27.9769265 ! 28Si
  real(kind=8),parameter :: s_mass= 31.9720712 ! 32S
  real(kind=8),parameter ::cl_mass= 34.9688527 ! 35Cl
  real(kind=8),parameter :: p_mass= 30.9737620 ! 31P
  real(kind=8),parameter ::br_mass= 78.9183376 ! 79Br
  real(kind=8),parameter :: i_mass=126.9044726 ! 127I
  real(kind=8),parameter :: b_mass= 11.0093052 ! 11B

  ! unit conversion factors
  real(kind=8),parameter :: auang = 0.5291772105 ! bohr and Angstrom
  real(kind=8),parameter :: aucm  =219474.631371 ! hartree and cm-1
  real(kind=8),parameter :: aukcal= 627.509474   ! hartree and kcal/mol
  real(kind=8),parameter :: auev  = 27.211386246 ! hartree and eV
  real(kind=8),parameter :: emass = 1822.88848   ! m_elec and amu

  ! other parameters
  real(kind=8),parameter :: pi = 3.14159265358979323844d0

end module constants
