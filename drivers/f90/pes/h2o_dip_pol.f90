!
!=============================================================================
!
! File: h2o_dip_pol.f90
! Original author: Tomislav Begusic
! Date: 03 April 2022
! Description: Simple test program to call stand-alone
!              Implements water dipole and polarizability 
!              functions for qTIP4P model following 
!              P. Hamm, J. Chem. Phys. 141, 184201 (2014).
!
! NOTES:
! (1) ATOMIC UNITS ARE USED THROUGHOUT
! (2) WE ASSUME THAT THE ATOMS ARE IN A LIST ARRANGED WITH MOLECULES
!     AS [ {O,H,H},{O,H,H},..... ].
!
! Disclaimer: please thoroughly check the output of this before 
! using it for anything interesting!
!
!=============================================================================
!


!
!====================================================================================
!
! The routines below implement dipole, its derivative, and polarizability of water. 
!
! Expects atoms in the order O H H O H H O H H ... 
!====================================================================================
!

! Dipole moment.
SUBROUTINE h2o_dipole(box, nat, atoms, dip, dip_der)

  IMPLICIT NONE
  
  DOUBLE PRECISION, INTENT(IN) :: box(3)
  INTEGER, INTENT(IN) :: nat
  DOUBLE PRECISION, INTENT(IN) :: atoms(nat, 3)
  DOUBLE PRECISION, INTENT(INOUT) :: dip(3)
  DOUBLE PRECISION, INTENT(INOUT) :: dip_der(nat, 3)

  DOUBLE PRECISION, PARAMETER :: pi = dacos(-1.0d0), twopi = 2 * pi
  DOUBLE COMPLEX, PARAMETER :: IU = CMPLX(0.0d0, 1.0d0, KIND=16)
  
  DOUBLE PRECISION, PARAMETER :: qm = -1.1128d0
  DOUBLE PRECISION, PARAMETER :: qh = - qm / 2
  DOUBLE PRECISION, PARAMETER :: charges(3) = (/ qm, qh, qh /) !(n_charged_atom_per_mol)
  DOUBLE PRECISION, PARAMETER :: gam = 0.73612d0, gam2 = 0.5d0 * (1-gam)
  DOUBLE PRECISION, PARAMETER :: a_iso = 1.47d0 * 1.88973d0**3
  DOUBLE PRECISION :: alpha(3) !(xyz)
  
  DOUBLE PRECISION :: dip_i(nat/3, 3) !Molecular dipole moments. (nmol, xyz)
  DOUBLE PRECISION :: dip_der_i(nat/3, nat, 3, 3) !Gradients of molecular dipole moments. (nmol, natoms, xyz, xyz)

  DOUBLE PRECISION :: atoms_charged(nat, 3) !List of atoms, with O replaced by M site ('charged atoms'). (natoms, xyz)
  DOUBLE PRECISION :: ro(nat/3, 3) !List of O atoms. (nmol, xyz)

  DOUBLE PRECISION :: E_stat(nat/3, 3) !(nmol, xyz)
  DOUBLE PRECISION :: T_tnsr(nat/3, nat/3, 3, 3, 3) !(nmol, nmol, n_charged_atom_per_mol, xyz, xyz)
  
  INTEGER :: i, iatom, k
  
  !Molecular polarizability.
  alpha = a_iso
  
  !Charged atoms (Site M instead of oxygen) and oxygen atoms stored separately.
  atoms_charged = atoms
  DO i=1, nat, 3
     !Apply PBC to hydrogens so that the molecule is not split over the box.
     !Probably not needed if we do not expect PBC applied to atoms before passed from i-pi.
     DO k = 1, 2
        atoms_charged(i + k, :) = atoms(i + k, :) - NINT((atoms(i + k, :) - atoms(i, :))/ box(:)) * box(:) 
     ENDDO
     !Compute M site position.
     atoms_charged(i, :) = gam * atoms(i, :) + gam2 * (atoms_charged(i+1, :) + atoms_charged(i+2, :))
     !Oxygen atoms.
     ro((i-1)/3 + 1, :) = atoms(i, :)
  ENDDO
  
  !Coordinate-dependent quantities: E_stat and T_tnsr. This is where we worry about Ewald summation.
  CALL calc_E_stat_and_T(E_stat, T_tnsr)

  DO i=1, nat/3
     iatom = 3*(i-1) + 1
     !Calculate dipole of i-th molecule. Using scaled charges for permanent dipole according to Hamm.
     dip_i(i, :) = calc_dipole(atoms_charged(iatom:iatom+2, :), charges / 1.3d0, E_stat(i, :))
     !Gradient of the i-th dipole moment.
     dip_der_i(i, :, : ,:) = calc_dipole_derivative(T_tnsr(i, :, :, :, :), i)
  ENDDO
  
  dip = SUM(dip_i, DIM=1)
  dip_der = SUM(dip_der_i(:, :, 3, :), DIM=1)

 CONTAINS

  SUBROUTINE calc_E_stat_and_T(E_stat, T_tnsr)

    DOUBLE PRECISION, INTENT(INOUT) :: E_stat(nat/3, 3), T_tnsr(nat/3, nat/3, 3, 3, 3)

    DOUBLE PRECISION :: a, rcut

    !----------------------------------------------
    ! Parameters for Ewald (calculated only once).
    !----------------------------------------------
    rcut = 1.5d0 * MINVAL(box) * MIN(0.5d0,1.2d0*nat**(-1.d0/6.d0))
    a = pi/rcut

    !----------------------------------------------
    ! Short-range part - sum over pairs.
    !----------------------------------------------
    CALL short_range_ew(atoms_charged, ro, a, rcut, E_stat, T_tnsr)

    !----------------------------------------------
    ! Long-range part - performs sum in k space.
    !----------------------------------------------
    CALL long_range_ew(atoms_charged, ro, a, E_stat)

  END SUBROUTINE calc_E_stat_and_T

  SUBROUTINE short_range_ew(r, ro, a, rcut, E_stat, T_tnsr)

    DOUBLE PRECISION, INTENT(IN) :: r(nat, 3), ro(nat/3, 3), a, rcut 
    DOUBLE PRECISION, INTENT(INOUT) :: E_stat(nat/3, 3), T_tnsr(nat/3, nat/3, 3, 3, 3)

    INTEGER :: i, j, k, jatom, nx, ny, nz, nmax
    DOUBLE PRECISION :: r_ij(3), dr, dr2, dr3, rcut2, screen, a_dr, gauss_part
    LOGICAL :: self_term

    nmax=NINT(rcut/MINVAL(box))
    E_stat = 0.0d0
    T_tnsr = 0.0d0
    rcut2 = rcut**2
    DO nx = -nmax, nmax
       DO ny = -nmax, nmax
          DO nz = -nmax, nmax
             DO i = 1, nat/3
                DO j = 1, nat/3
                   DO k = 1, 3
                      jatom = 3 * (j - 1) + k
                      r_ij(:) = ro(i,:) - r(jatom, :)
                      r_ij(:) = r_ij(:) - box(:) * (NINT(r_ij(:)/box(:)) + (/nx, ny, nz/) )
                      dr2 = SUM(r_ij**2)
                      IF (dr2 .LT. rcut2) THEN
                         dr = SQRT(dr2)
                         dr3 = dr*dr2
                         self_term = (i .EQ. j) .AND. (nx .EQ. 0) .AND. (ny .EQ. 0).AND. (nz .EQ. 0)
                         a_dr = a * dr
                         gauss_part = 2.0d0 / SQRT(pi) * a_dr * EXP(-a_dr**2)
                         screen = short_range_ew_screen(a_dr, gauss_part, self_term)
                         E_stat(i, :) = E_stat(i, :) + charges(k) * r_ij(:) / dr3 * screen
                         T_tnsr(i, j, k, :, :) = T_tnsr(i, j, k, :, :) + charges(k) * short_range_T_tnsr(r_ij, dr2, dr3, a_dr, gauss_part, screen)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE short_range_ew

  FUNCTION short_range_ew_screen(a_r, gauss_part, self_term) RESULT(sc)

    DOUBLE PRECISION :: sc

    DOUBLE PRECISION, INTENT(IN) :: a_r, gauss_part
    LOGICAL, INTENT(IN) :: self_term

    sc = erfc(a_r) + gauss_part
    IF (self_term) sc = sc - 1.0d0 !Self-interaction term.

  END FUNCTION short_range_ew_screen

  FUNCTION short_range_T_tnsr(r_ij, r2, r3, a_r, gauss_part, screen) RESULT(T_ij)

    DOUBLE PRECISION :: T_ij(3, 3)

    DOUBLE PRECISION, INTENT(IN) :: r_ij(3), r2, r3, a_r, gauss_part, screen

    INTEGER :: l

    T_ij = - outer(r_ij, r_ij) / (r3 * r2) * (3.0d0 * screen + 2.0d0 * a_r**2 * gauss_part)
    DO l = 1, 3
       T_ij(l, l) = T_ij(l, l) + screen / r3
    ENDDO

  END FUNCTION short_range_T_tnsr

  SUBROUTINE long_range_ew(r, ro, a, E_stat)

    DOUBLE PRECISION, INTENT(IN) :: r(nat, 3), ro(nat/3, 3), a    
    DOUBLE PRECISION, INTENT(INOUT) :: E_stat(nat/3, 3)

    INTEGER :: l, kx, ky, kz, kmax
    DOUBLE PRECISION :: b, f, rk(3), rk2, rkmax2, lat(3), q(nat), tmp(nat/3)
    DOUBLE COMPLEX :: sk

    kmax = INT(a * MAXVAL(box))
    lat(:) = twopi/box(:) 
    rkmax2 = (twopi * a)**2
    b = 0.25d0/(a**2)
    f = 4.0d0 * pi / PRODUCT(box)     ! 4pi / V !
    q(:) = PACK(SPREAD(charges(:), 2, nat/3), .TRUE.)

    DO kx = 0,kmax
       IF (kx .EQ. 1) f = 2.d0*f
       rk(1) = lat(1)*kx
       DO ky = -kmax,kmax
          rk(2) = lat(2)*ky
          DO kz = -kmax,kmax
             rk(3) = lat(3)*kz
             rk2 = SUM(rk(:)**2)
             IF (rk2 .LT. rkmax2 .AND. rk2 .GT. EPSILON(0.d0)) THEN
                sk = f * (EXP(-b*rk2) / rk2) * SUM(q * EXP(-IU * k_dot_r(nat, rk, r)))
                tmp = AIMAG(EXP(IU * k_dot_r(nat/3, rk, ro)) * sk)
                DO l = 1, 3
                   E_stat(:,l) = E_stat(:, l) + rk(l)  * tmp(:)
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE long_range_ew

  FUNCTION calc_dipole(atoms, charges, E_stat)

    DOUBLE PRECISION, INTENT(IN) :: atoms(3,3), charges(3), E_stat(3)

    DOUBLE PRECISION :: calc_dipole(3)

    INTEGER :: k

    calc_dipole = 0.0d0
    !Permanent dipoles.
    DO k = 1, 3
       calc_dipole(:) = calc_dipole(:) + charges(k) * atoms(k, :) 
    ENDDO
    !Induced dipoles from electrostatic interaction with charges on other water molecules.
    calc_dipole(:) = calc_dipole(:) + alpha * E_stat

  END FUNCTION calc_dipole

  FUNCTION calc_dipole_derivative(T_tnsr, i) RESULT(dip_der)

    DOUBLE PRECISION, INTENT(IN) :: T_tnsr(nat/3, 3, 3, 3)
    INTEGER :: i

    DOUBLE PRECISION :: dip_der(nat, 3, 3)

    INTEGER :: iatom, j, jatom, k

    dip_der = 0.0d0
    iatom = 3 * (i - 1) + 1
    DO j = 1, nat/3
       jatom = 3*(j-1) + 1
       !Sum for i-th oxygen atom.
       DO k = 1, 3
          dip_der(iatom, :, :) =  dip_der(iatom, :, :) + T_tnsr(j, k, :, :)
       ENDDO
       !j-th molecule oxygen atom. Includes i==j term.
       dip_der(jatom, :, :) = dip_der(jatom, :, :) - gam * T_tnsr(j, 1, :, :)
       !j-th molecule hydrogen atoms.
       DO k = 2, 3
          dip_der(jatom + k - 1, :, :) = - gam2 * T_tnsr(j, 1, :, :) - T_tnsr(j, k, :, :)
       ENDDO
    ENDDO
    DO k = 1, 3
       !Multiply induced part of dipole derivative by alpha.
       DO j = 1, nat
          dip_der(j, k, :) = alpha(k) * dip_der(j, k, :)
       ENDDO
       !Add gradient of the permanent molecular dipole moment.
       dip_der(iatom, k, k) = dip_der(iatom, k, k)   + gam * qm !i-th oxygen.
       dip_der(iatom+1:iatom+2, k, k) = qh + gam2 * qm !i-th hydrogen.
    ENDDO

  END FUNCTION calc_dipole_derivative


  FUNCTION k_dot_r(n, k, r)

    DOUBLE PRECISION :: k_dot_r(n)

    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: k(3), r(n, 3)

    INTEGER :: i

    k_dot_r = 0.0d0
    DO i = 1, 3
       k_dot_r(:) = k_dot_r(:) + k(i) * r(:, i) 
    ENDDO

  END FUNCTION k_dot_r

  FUNCTION outer(a, b)

    DOUBLE PRECISION :: a(3), b(3)

    DOUBLE PRECISION :: outer(3,3)

    INTEGER :: i

    DO i=1, 3
       outer(i, :)  = a(i) * b(:)
    ENDDO

  END FUNCTION outer

END SUBROUTINE h2o_dipole

