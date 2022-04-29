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
SUBROUTINE h2o_dipole(box, nat, atoms, compute_der, dip, dip_der, pol)

  IMPLICIT NONE
  
  DOUBLE PRECISION, INTENT(IN) :: box(3)
  INTEGER, INTENT(IN) :: nat
  DOUBLE PRECISION, INTENT(IN) :: atoms(nat, 3)
  LOGICAL, INTENT(IN) :: compute_der
  DOUBLE PRECISION, INTENT(INOUT) :: dip(3)
  DOUBLE PRECISION, INTENT(INOUT) :: dip_der(nat, 3)
  DOUBLE PRECISION, INTENT(INOUT) :: pol(3, 3)

  DOUBLE PRECISION, PARAMETER :: pi = DACOS(-1.0d0), twopi = 2 * pi, sqrtpi = SQRT(pi)
  DOUBLE COMPLEX, PARAMETER :: IU = CMPLX(0.0d0, 1.0d0, KIND=16)
  
  DOUBLE PRECISION, PARAMETER :: qm = -1.1128d0
  DOUBLE PRECISION, PARAMETER :: qh = - qm / 2
  DOUBLE PRECISION, PARAMETER :: charges(3) = (/ qm, qh, qh /) !(n_charged_atom_per_mol)
  DOUBLE PRECISION, PARAMETER :: gam = 0.73612d0, gam2 = 0.5d0 * (1-gam)
  DOUBLE PRECISION, PARAMETER :: a_iso = 1.47d0 * 1.88973d0**3

  DOUBLE PRECISION :: alpha(3) !(xyz)
  DOUBLE PRECISION :: atoms_charged(nat, 3) !List of atoms, with O replaced by M site ('charged atoms'). (natoms, xyz)
  DOUBLE PRECISION :: ro(nat/3, 3) !List of O atoms. (nmol, xyz)
  DOUBLE PRECISION :: dip_der_full(nat, 3, 3) !Gradient of dipole moment. (natoms, xyz, xyz)

  DOUBLE PRECISION :: dip_ind(3) !(nmol, xyz)
  DOUBLE PRECISION :: dip_ind_der(nat, 3, 3) !(nmol, nat, xyz, xyz)
  DOUBLE PRECISION :: pol_ind(3, 3) !(nmol, nmol, xyz, xyz)
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
  
  !Coordinate-dependent quantities: dip_ind and dip_ind_der. This is where we worry about Ewald summation.
  CALL calc_induced_part(dip_ind, dip_ind_der, pol_ind)

  !Calculate total dipole moment. Using scaled charges for permanent dipole according to Hamm.
  dip(:) = calc_dipole(atoms_charged, charges / 1.3d0, dip_ind)
  IF (compute_der) THEN
     dip_der_full(:, : ,:) = calc_dipole_derivative(charges / 1.3d0, dip_ind_der)
  ELSE
     dip_der_full(:, :, :) = 0.0d0
  ENDIF
  dip_der = dip_der_full(:, 3, :)
  pol = calc_polarizability(pol_ind)

 CONTAINS

  SUBROUTINE calc_induced_part(dip_ind, dip_ind_der, pol_ind)

    DOUBLE PRECISION, INTENT(INOUT) :: dip_ind(3), dip_ind_der(nat, 3, 3), pol_ind(3, 3)

    DOUBLE PRECISION :: a, rcut

    !----------------------------------------------
    ! Parameters for Ewald (calculated only once).
    !----------------------------------------------
    rcut = 1.5d0 * MINVAL(box) * MIN(0.5d0,1.2d0*nat**(-1.d0/6.d0))
    a = pi/rcut

    !----------------------------------------------
    ! Short-range part - sum over pairs.
    !----------------------------------------------
    CALL short_range_ew(atoms_charged, ro, a, rcut, dip_ind, dip_ind_der, pol_ind)

    !----------------------------------------------
    ! Long-range part - performs sum in k space.
    !----------------------------------------------
    CALL long_range_ew(atoms_charged, ro, a, dip_ind, dip_ind_der, pol_ind)

  END SUBROUTINE calc_induced_part

  SUBROUTINE short_range_ew(r, ro, a, rcut, dip_ind, dip_ind_der, pol_ind)

    DOUBLE PRECISION, INTENT(IN) :: r(nat, 3), ro(nat/3, 3), a, rcut 
    DOUBLE PRECISION, INTENT(INOUT) :: dip_ind(3), dip_ind_der(nat, 3, 3), pol_ind(3, 3)

    INTEGER :: i, j, k, l, jatom, nx, ny, nz, nmax
    DOUBLE PRECISION :: r_ij_0(3), r_ij(3), dr, dr2, dr3, rcut2, screen, a_dr, gauss_part, T_tnsr(3, 3)
    LOGICAL :: self_term

    nmax=NINT(rcut/MINVAL(box))
    dip_ind = 0.0d0
    dip_ind_der = 0.0d0
    rcut2 = rcut**2
    DO i = 1, nat/3
       iatom = 3 * (i - 1) + 1
       DO j = 1, nat/3
          jatom = 3 * (j - 1)
          DO k = 1, 3
             T_tnsr = 0.0d0
             !Nearest neighbor.
             r_ij_0(:) = ro(i,:) - r(jatom + k, :)
             r_ij_0(:) = r_ij_0(:) - box(:) * NINT(r_ij_0(:)/box(:))
             DO nx = -nmax, nmax
                DO ny = -nmax, nmax
                   DO nz = -nmax, nmax
                      r_ij(:) = r_ij_0(:) - box(:) * (/nx, ny, nz/) 
                      dr2 = SUM(r_ij**2)
                      IF (dr2 .LT. rcut2) THEN
                         ! Helper variables that are used more than once or convenient to compute separately.
                         dr = SQRT(dr2)
                         dr3 = dr*dr2
                         self_term = (i .EQ. j) .AND. (nx .EQ. 0) .AND. (ny .EQ. 0).AND. (nz .EQ. 0)
                         a_dr = a * dr
                         gauss_part = 2.0d0 / sqrtpi * a_dr * EXP(-a_dr**2)
                         screen = short_range_ew_screen(a_dr, gauss_part, self_term)
                         ! Contribution to the T tensor for given i, j, k, n (used for dipole derivative).
                         IF (compute_der) THEN
                            T_tnsr = T_tnsr + short_range_T_tnsr(r_ij, dr2, dr3, a_dr, gauss_part, screen)
                         ENDIF
                         !--------------------------------------------------------------------------------------
                         ! Contribution to the electric field of molecule i from atom k in molecule j in cell n.
                         !--------------------------------------------------------------------------------------
                         dip_ind(:) = dip_ind(:) + alpha * charges(k) * r_ij(:) / dr3 * screen
                         !--------------------------------------------------------------------------------------
                         !--------------------------------------------------------------------------------------
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             !--------------------------------------------------------------------------------------------------------
             ! Compute the derivative of the electric field of molecule i w.r.t. coordinates of atom k of molecule j.
             !--------------------------------------------------------------------------------------------------------
             IF (compute_der) THEN
                DO l = 1, 3
                   T_tnsr(l, :) =  alpha(l) * T_tnsr(l, :)
                ENDDO
                T_tnsr =  charges(k) * T_tnsr
                !Sum for i-th oxygen atom.
                dip_ind_der(iatom, :, :) =  dip_ind_der(iatom, :, :) + T_tnsr(:, :)
                !Derivative of electric field at molecule i w.r.t. coordinates of molecule j (includes self-term i==j).
                CALL dip_ind_der_ij(dip_ind_der(jatom + 1: jatom + 3, :, :), T_tnsr(:, :), k)
             ENDIF
             !--------------------------------------------------------------------------------------------------------
             !--------------------------------------------------------------------------------------------------------
          ENDDO
       ENDDO
    ENDDO

    
    !rcut2 = 9 * rcut2
    !nmax=NINT(3*rcut/MINVAL(box))
    pol_ind = 0.0d0
    DO i = 1, nat/3
       DO j = 1, nat/3
          r_ij_0(:) = ro(i,:) - ro(j, :)
          r_ij_0(:) = r_ij_0(:) - box(:) * NINT(r_ij_0(:)/box(:))
          T_tnsr = 0.0d0
          DO nx = -nmax, nmax
             DO ny = -nmax, nmax
                DO nz = -nmax, nmax
                   r_ij(:) = r_ij_0(:) - box(:) * (/nx, ny, nz/) 
                   dr2 = SUM(r_ij**2)
                   self_term = (i .EQ. j) .AND. (nx .EQ. 0) .AND. (ny .EQ. 0).AND. (nz .EQ. 0)
                   IF (dr2 .LT. rcut2 .AND. (.NOT. self_term)) THEN
                      dr = SQRT(dr2)
                      dr3 = dr*dr2
                      T_tnsr = T_tnsr + short_range_T_tnsr(r_ij, dr2, dr3, 0.0d0, 0.0d0, 1.0d0)
                      !a_dr = a * dr
                      !gauss_part = 2.0d0 / sqrtpi * a_dr * EXP(-a_dr**2)
                      !screen = short_range_ew_screen(a_dr, gauss_part, .FALSE.)
                      !T_tnsr = T_tnsr + short_range_T_tnsr(r_ij, dr2, dr3, a_dr, gauss_part, screen)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          !IF (i .EQ. j) THEN
          !   DO k = 1, 3
          !      T_tnsr(k, k) = T_tnsr(k, k) - 4.0d0 * a**3 / (3.0d0 * sqrtpi)
          !   ENDDO
          !ENDIF
          DO k = 1, 3
             pol_ind(k, :) = pol_ind(k, :) + alpha(k) * T_tnsr(k, :) * alpha(:)
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

  SUBROUTINE long_range_ew(r, ro, a, dip_ind, dip_ind_der, pol_ind)

    DOUBLE PRECISION, INTENT(IN) :: r(nat, 3), ro(nat/3, 3), a    
    DOUBLE PRECISION, INTENT(INOUT) :: dip_ind(3), dip_ind_der(nat, 3, 3), pol_ind(3, 3)

    INTEGER :: i, k, l, kx, ky, kz, kmax
    DOUBLE PRECISION :: b, f, rk(3), rk_out(3, 3), rk2, rkmax2, lat(3), q(nat), tnsr_tmp(3, 3)
    DOUBLE COMPLEX :: sk_i(nat), exp_ikr(nat/3), sk, sk_o(3)

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
                sk_i = f * (EXP(-b*rk2) / rk2) * q * EXP(-IU * k_dot_r(nat, rk, r))
                exp_ikr = EXP(IU * k_dot_r(nat/3, rk, ro))
                sk = SUM(sk_i)
                sk_o = alpha * SUM(exp_ikr)
                dip_ind(:) = dip_ind(:) + rk(:)  * AIMAG(sk_o * sk)
                IF (compute_der) THEN
                   rk_out = outer(rk, rk)
                   DO i = 1, nat/3
                      iatom = 3 * (i-1) + 1
                      !Sum for i-th oxygen atom.
                      DO l = 1, 3
                         dip_ind_der(iatom, l, :) = dip_ind_der(iatom, l, :) + alpha(l) * rk_out(l, :) * REAL(exp_ikr(i) * sk, KIND=KIND(1.0d0))
                      ENDDO
                      !Derivatives with respect to atom k of molecule j.
                      DO k = 1, 3
                         DO l = 1, 3
                            tnsr_tmp(l, :) = REAL(sk_o(l) * sk_i(iatom + k - 1), KIND=KIND(1.0d0)) * rk_out(l, :)
                         ENDDO
                         CALL dip_ind_der_ij(dip_ind_der(iatom: iatom + 2, :, :), tnsr_tmp, k)
                      ENDDO
                   ENDDO
                ENDIF
                !DO l = 1, 3
                !   pol_ind(l, :) = pol_ind(l, :) + f * EXP(-b*rk2) / rk2 * ABS(sk_o(l)) * rk_out(l, :) * ABS(sk_o(:))
                !ENDDO
             ENDIF
             !IF (rk2 .LT. EPSILON(0.0d0)) THEN
             !   DO l = 1, 3
             !      pol_ind(l, l) = pol_ind(l, l) + 4.0d0 * pi / PRODUCT(box) * (nat / 3.0d0 * alpha(l))**2
             !   ENDDO
             !ENDIF 
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE long_range_ew

  SUBROUTINE dip_ind_der_ij(dEi_drj, tnsr, k)

    DOUBLE PRECISION, INTENT(INOUT) :: dEi_drj(3, 3, 3)
    DOUBLE PRECISION, INTENT(IN) :: tnsr(3, 3)
    INTEGER, INTENT(IN) :: k

    INTEGER :: h

    IF (k .EQ. 1) THEN
       !j-th molecule oxygen atom. Includes i==j term.
       dEi_drj(k, :, :) = dEi_drj(k, :, :) - gam * tnsr(:, :)
       !j-th molecule hydrogen atoms (part that comes from M site).
       DO h = 1, 2
          dEi_drj(k + h, :, :) = dEi_drj(k + h, :, :) - gam2 * tnsr(:, :)
       ENDDO
    ELSE
       !j-th molecule hydrogen atoms (part that comes from hydrogen atom positions).
       dEi_drj(k, :, :) = dEi_drj(k, :, :) - tnsr(:, :)
    ENDIF   

  END SUBROUTINE dip_ind_der_ij

  FUNCTION calc_dipole(atoms, charges, dip_ind)

    DOUBLE PRECISION, INTENT(IN) :: atoms(nat,3), charges(3), dip_ind(3)

    DOUBLE PRECISION :: calc_dipole(3)

    INTEGER :: i, iatom, k

    calc_dipole = 0.0d0
    !Permanent dipoles.
    DO i = 1, nat/3
       iatom = 3 * (i - 1)
       DO k = 1, 3
          calc_dipole(:) = calc_dipole(:) + charges(k) * atoms(iatom + k, :) 
       ENDDO
       !Induced dipoles from electrostatic interaction with charges on other water molecules.
    ENDDO
    calc_dipole(:) = calc_dipole(:) + dip_ind(:)

  END FUNCTION calc_dipole

  FUNCTION calc_dipole_derivative(charges, dip_ind_der) RESULT(dip_der)

    DOUBLE PRECISION, INTENT(IN) :: charges(3)
    DOUBLE PRECISION, INTENT(IN) :: dip_ind_der(nat, 3, 3)

    DOUBLE PRECISION :: dip_der(nat, 3, 3)

    INTEGER :: iatom, i, k

    dip_der = 0.0d0
    DO i=1, nat/3
       iatom = 3 * (i - 1) + 1
       DO k = 1, 3
          !Gradient of the permanent molecular dipole moment.
          dip_der(iatom, k, k) = dip_der(iatom, k, k) + gam * charges(1) !i-th oxygen.
          dip_der(iatom+1:iatom+2, k, k) = dip_der(iatom+1:iatom+2, k, k) + charges(2) + gam2 * charges(1) !i-th hydrogen.
       ENDDO
    ENDDO
    !Multiply induced part of dipole derivative by alpha.
    DO k = 1, 3
       dip_der(:, k, :) = dip_der(:, k, :) + dip_ind_der(:, k, :)
    ENDDO

  END FUNCTION calc_dipole_derivative

  FUNCTION calc_polarizability(pol_ind) RESULT(pol)

    DOUBLE PRECISION, INTENT(IN) :: pol_ind(3,3)

    DOUBLE PRECISION :: pol(3, 3)

    INTEGER :: k

    pol = 0.0d0
    DO k = 1, 3
    !   pol(k, k) = pol(k, k) + alpha(k)
    ENDDO
    pol = pol + pol_ind 

  END FUNCTION calc_polarizability

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

