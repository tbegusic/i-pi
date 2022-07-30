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
  DOUBLE PRECISION, PARAMETER :: a_iso = 1.47d0 * 1.88973d0**3 !Isotropic polarizability from Hamm's paper.
  !Anisotropic polarizability from Hamm's paper.
  DOUBLE PRECISION, PARAMETER :: a_aniso(3) = (/ 1.626d0, 1.495d0, 1.286d0 /) * 1.88973d0**3
  !Anisotropic polarizability from J. Comput. Chem. 2016, 37, 2125–2132.
  !Was not better than the one from Hamm.
  !DOUBLE PRECISION, PARAMETER :: a_aniso(3) = (/ 1.37071d0, 1.41205d0, 1.46317d0 /) * 1.88973d0**3

  INTEGER :: nmol

  DOUBLE PRECISION :: alpha(nat/3, 3, 3) !(nmol, xyz, xyz)
  DOUBLE PRECISION :: atoms_charged(nat, 3) !List of atoms, with O replaced by M site ('charged atoms'). (natoms, xyz)
  DOUBLE PRECISION :: ro(nat/3, 3) !List of O atoms. (nmol, xyz)
  DOUBLE PRECISION :: dip_der_full(nat, 3, 3) !Gradient of dipole moment. (natoms, xyz, xyz)
  DOUBLE PRECISION :: dalpha_dr(nat, 3, 3, 3) !Gradient of molecular polarizabilities. (natoms, xyz, xyz, xyz)

  DOUBLE PRECISION :: dip_ind(3) !(nmol, xyz)
  DOUBLE PRECISION :: dip_ind_der(nat, 3, 3) !(nmol, nat, xyz, xyz)
  DOUBLE PRECISION :: pol_ind(3, 3) !(nmol, nmol, xyz, xyz)
  INTEGER :: i, k
  
  !Number of molecules.
  nmol = nat / 3

  !Molecular polarizability.
  CALL calc_alpha
  
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

  SUBROUTINE calc_alpha

    INTEGER :: i, iatom, j, k
    DOUBLE PRECISION :: O(3, 3), dO_dr(3, 3, 3, 3), norm_tmp, tnsr_tmp(3, 3), dist_vec(3)

    alpha = 0.0d0
    DO i = 1, nmol

       iatom = 3 * (i - 1) + 1
       DO k = 1, 3
          !alpha(i, k, k) = a_iso              !Isotropic molecular polarizability.
          alpha(i, k, k) = a_aniso(k)          !Anisotropic molecular polarizability.
       ENDDO
       !alpha(i, :, :) = calc_aniso_alpha(i, atoms(iatom:iatom+2, :)) !Anisotropic molecular polarizability.

       !x-axis: (rH1 - rH2) normalized.
       dist_vec = atoms(iatom + 1, :) - atoms(iatom + 2, :)
       norm_tmp = NORM2(dist_vec)
       O(:, 1) = dist_vec / norm_tmp
       IF (compute_der) THEN
          tnsr_tmp = rot_grad_tnsr(O(:, 1), norm_tmp)
          dO_dr(1, :, :, 1) = 0.0d0
          dO_dr(2, :, :, 1) = tnsr_tmp
          dO_dr(3, :, :, 1) = -tnsr_tmp
       ENDIF
       !y-axis: rO - rH1 is orthonormalized w.r.t. rx.
       dist_vec = atoms(iatom, :) - atoms(iatom + 1, :)
       O(:, 2) = dist_vec - DOT_PRODUCT(dist_vec, O(:, 1)) * O(:, 1)
       norm_tmp = NORM2(O(:, 2))
       O(:, 2) = O(:, 2) / norm_tmp
       IF (compute_der) THEN
          tnsr_tmp = rot_grad_tnsr(O(:, 2), norm_tmp)
          dO_dr(1, :, :, 2) = MATMUL(rot_grad_tnsr(O(:, 1), 1.0d0), tnsr_tmp)
          dO_dr(3, :, :, 2) = - MATMUL(outer(MATMUL(dO_dr(3, :, :, 1), dist_vec), O(:, 1)) + DOT_PRODUCT(O(:, 1), dist_vec) * dO_dr(3, :, :, 1), tnsr_tmp)
          dO_dr(2, :, :, 2) = -dO_dr(3, :, :, 2) - dO_dr(1, :, :, 2)
       ENDIF
       !z-axis: Cross product between rx and ry (perpendicular to the molecular plane).
       O(:, 3) = cross_product(O(:, 1), O(:, 2))
       IF (compute_der) THEN
          DO j = 1, 3
             DO k = 1, 3
                dO_dr(j, k, :, 3) = cross_product(dO_dr(j, k, :, 1), O(:, 2)) + cross_product(O(:, 1), dO_dr(j, k, :, 2)) 
             ENDDO
          ENDDO
       ENDIF

       !Rotated polarizability and its gradient.
       IF (compute_der) THEN
          DO j = 1, 3
             DO k = 1, 3
                tnsr_tmp = MATMUL(MATMUL(dO_dr(j, k, :, :), alpha(i, :, :)), TRANSPOSE(O))
                dalpha_dr(iatom + j - 1, k, :, :) = tnsr_tmp + TRANSPOSE(tnsr_tmp)
             ENDDO
          ENDDO
       ENDIF
       alpha(i, :, :) = MATMUL(MATMUL(O, alpha(i, :, :)), TRANSPOSE(O))
    ENDDO

  END SUBROUTINE calc_alpha

  ! Computes matrix (Id - outer(vec, vec)) / norm, which appears in the
  ! computation of the gradient of the rotation matrix.
  FUNCTION rot_grad_tnsr(vec, norm) RESULT(tnsr)

    DOUBLE PRECISION, INTENT(IN) :: vec(3), norm

    DOUBLE PRECISION :: tnsr(3, 3)

    INTEGER :: k

    tnsr = - outer(vec, vec)
    DO k = 1, 3
       tnsr(k, k) = tnsr(k, k) + 1.0d0
    ENDDO
    tnsr = tnsr / norm

  END FUNCTION rot_grad_tnsr

  SUBROUTINE calc_induced_part(dip_ind, dip_ind_der, pol_ind)

    DOUBLE PRECISION, INTENT(INOUT) :: dip_ind(3), dip_ind_der(nat, 3, 3), pol_ind(3, 3)

    DOUBLE PRECISION :: a, rcut

    !----------------------------------------------
    ! Parameters for Ewald (calculated only once).
    !----------------------------------------------
    rcut = 1.3d0 * MINVAL(box) * MIN(0.5d0,1.2d0*nat**(-1.d0/6.d0))
    a = 1.3d0 * pi/rcut

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

    DOUBLE PRECISION, INTENT(IN) :: r(nat, 3), ro(nmol, 3), a, rcut 
    DOUBLE PRECISION, INTENT(INOUT) :: dip_ind(3), dip_ind_der(nat, 3, 3), pol_ind(3, 3)

    INTEGER :: i, j, k, l, m, iatom, jatom, latom, nx, ny, nz, nmax
    DOUBLE PRECISION :: r_ij_0(3), r_ij(3), dr, dr2, dr3, rcut2, screen, a_dr, gauss_part
    DOUBLE PRECISION :: T_tnsr(3, 3), T_tnsr_sum(3, 3), E_ij(3), a3_self_term, prefac
    LOGICAL :: self_term

    nmax=NINT(rcut/MINVAL(box))
    dip_ind = 0.0d0
    dip_ind_der = 0.0d0
    rcut2 = rcut**2
    prefac = 2.0d0 / sqrtpi
    DO i = 1, nmol
       iatom = 3 * (i - 1) + 1
       T_tnsr_sum = 0.0d0
       E_ij = 0.0d0
       DO j = 1, nmol
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
                         self_term = (i .EQ. j) .AND. (nx .EQ. 0) .AND. (ny .EQ. 0).AND. (nz .EQ. 0)
                         CALL compute_helper_vars(a, dr2, prefac, self_term, dr, dr3, a_dr, gauss_part, screen)
                         ! Contribution to the electric field on molecule i induced by atom k in molecule j of cell n.
                         E_ij(:) = E_ij(:) + charges(k) * r_ij(:) / dr3 * screen
                         ! Contribution to the T tensor for given i, j, k, n (used for dipole derivative).
                         IF (compute_der) T_tnsr = T_tnsr + short_range_T_tnsr(r_ij, dr2, dr3, a_dr, gauss_part, screen)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             !--------------------------------------------------------------------------------------------------------
             ! Compute the derivative of the electric field of molecule i w.r.t. coordinates of atom k of molecule j.
             !--------------------------------------------------------------------------------------------------------
             IF (compute_der) THEN
                T_tnsr =  charges(k) * MATMUL(alpha(i, :, :), T_tnsr)
                !Derivative of induced dipole at molecule i w.r.t. coordinates of molecule j (includes self-term i==j).
                CALL dip_ind_der_ij(dip_ind_der(jatom + 1: jatom + 3, :, :), T_tnsr(:, :), k)
                !Sum for derivative w.r.t. i-th molecule atoms.
                T_tnsr_sum = T_tnsr_sum + T_tnsr
             ENDIF
             !--------------------------------------------------------------------------------------------------------
             !--------------------------------------------------------------------------------------------------------
          ENDDO
       ENDDO
       !----------------------------------------------------------------------------------------------------------
       ! Contribution to the dipole moment of molecule i induced by atom k in molecule j summed over all cells.
       !----------------------------------------------------------------------------------------------------------
       dip_ind(:) = dip_ind(:) + MATMUL(alpha(i, :, :), E_ij(:))
       !----------------------------------------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------

       !--------------------------------------------------------------------------------------------------------
       ! Compute the derivative of the electric field of molecule i w.r.t. coordinates of atoms in molecule i.
       !--------------------------------------------------------------------------------------------------------
       IF (compute_der) THEN
          !Oxygen of molecule i.
          dip_ind_der(iatom, :, :) = dip_ind_der(iatom, :, :) + T_tnsr_sum(:, :)
          !Derivative of alpha w.r.t. atoms of i-th molecule (loop over oxygen and hydrogens, and over xyz).
          DO l = 1, 3
             latom = iatom + l - 1
             DO m = 1, 3
                dip_ind_der(latom, m, :) = dip_ind_der(latom, m, :) + MATMUL(dalpha_dr(latom, :, m, :), E_ij)
             ENDDO
          ENDDO
       ENDIF
       !--------------------------------------------------------------------------------------------------------
       !--------------------------------------------------------------------------------------------------------
    ENDDO

    !--------------------------------------------------------------------------------------------------------
    ! Compute induced part of the polarizability.
    !--------------------------------------------------------------------------------------------------------
    a3_self_term = 4.0d0 * a**3 / (3.0d0 * sqrtpi)
    pol_ind = 0.0d0
    DO i = 1, nmol
       DO j = i, nmol
          r_ij_0(:) = ro(i,:) - ro(j, :)
          r_ij_0(:) = r_ij_0(:) - box(:) * NINT(r_ij_0(:)/box(:))
          T_tnsr = 0.0d0
          !Sum over cells and consider only non-self terms.
          DO nx = -nmax, nmax
             DO ny = -nmax, nmax
                DO nz = -nmax, nmax
                   r_ij(:) = r_ij_0(:) - box(:) * (/nx, ny, nz/) 
                   dr2 = SUM(r_ij**2)
                   self_term = (i .EQ. j) .AND. (nx .EQ. 0) .AND. (ny .EQ. 0).AND. (nz .EQ. 0)
                   IF (dr2 .LT. rcut2 .AND. (.NOT. self_term)) THEN
                      CALL compute_helper_vars(a, dr2, prefac, .FALSE., dr, dr3, a_dr, gauss_part, screen)
                      T_tnsr = T_tnsr + short_range_T_tnsr(r_ij, dr2, dr3, a_dr, gauss_part, screen)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          !Multiply computed tensor by alpha and add to system's polarizability.
          IF (i .EQ. j) THEN
             !Self term, include one for each value of i. 
             pol_ind = pol_ind - a3_self_term * MATMUL(alpha(i, :, :), alpha(j, :, :))
          ELSE
             !Not a self term, so we have to add also the transpose, to correct for including only j>i.
             ! alpha_i * T_ij * alpha_j = (alpha_j * T_ji * alpha_i)^T
             T_tnsr = MATMUL(MATMUL(alpha(i, :, :), T_tnsr), alpha(j, :, :))
             pol_ind = pol_ind + T_tnsr + TRANSPOSE(T_tnsr)
          ENDIF
       ENDDO
    ENDDO
    !--------------------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------

  END SUBROUTINE short_range_ew

  SUBROUTINE compute_helper_vars(a, dr2, prefac, self_term, dr, dr3, a_dr, gauss_part, screen)

    DOUBLE PRECISION, INTENT(IN) :: a, dr2, prefac
    LOGICAL, INTENT(IN) :: self_term
    DOUBLE PRECISION, INTENT(OUT) :: dr, dr3, a_dr, gauss_part, screen

    ! Helper variables for real-space summation:
    dr = SQRT(dr2)                                                  !Absolute value of distance.
    dr3 = dr*dr2                                                    !dr cubed.
    a_dr = a * dr                                                   !a * dr.
    gauss_part = prefac * a_dr * EXP(-a_dr**2)                      !Gaussian multiplied by a_dr and the prefactor.
    screen = short_range_ew_screen(a_dr, gauss_part, self_term)     !screen term equal to erfc(x) + gauss_part (-1 if self term).

  END SUBROUTINE compute_helper_vars

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

    DOUBLE PRECISION, INTENT(IN) :: r(nat, 3), ro(nmol, 3), a    
    DOUBLE PRECISION, INTENT(INOUT) :: dip_ind(3), dip_ind_der(nat, 3, 3), pol_ind(3, 3)

    INTEGER :: i, l, m, iatom, latom, kx, ky, kz, kmax
    DOUBLE PRECISION :: b, f, f_mod, rk(3), rk_out(3, 3), rk2, rkmax2, lat(3), q(nat), tnsr_tmp(3, 3), prefac, sum_alpha(3, 3) 
    DOUBLE COMPLEX :: sk_i(nat), exp_ikr(nmol), sk, sk_o(3, 3), sk_o_times_rk(3)

    kmax = INT(1.3d0*a * MAXVAL(box))
    lat(:) = twopi/box(:) 
    rkmax2 = (1.3d0*twopi * a)**2
    b = 0.25d0/(a**2)
    f = 4.0d0 * pi / PRODUCT(box) ! 4pi / V !
    f_mod = f
    q(:) = PACK(SPREAD(charges(:), 2, nmol), .TRUE.)

    DO kx = 0,kmax
       IF (kx .EQ. 1) f_mod = 2.d0*f
       rk(1) = lat(1)*kx
       DO ky = -kmax,kmax
          rk(2) = lat(2)*ky
          DO kz = -kmax,kmax
             rk(3) = lat(3)*kz
             rk2 = SUM(rk(:)**2)
             IF (rk2 .LT. rkmax2 .AND. rk2 .GT. EPSILON(0.d0)) THEN
                !Helper variables:
                prefac = f_mod * EXP(-b*rk2) / rk2         !Prefactor common to all properties.
                sk_i = q * EXP(-IU * k_dot_r(nat, rk, r))  !charge(i) * exp(i k r_i), where r_i are charged sites.
                sk = SUM(sk_i)                             !Structure factor for the charged sites.
                exp_ikr = EXP(IU * k_dot_r(nmol, rk, ro))  !exp(i k r_o), where r_o are oxygen atoms.
                sk_o = 0.0d0                               !Structure factor for the oxygen atoms.
                DO i = 1, nmol
                   sk_o = sk_o + alpha(i, :, :) * exp_ikr(i) 
                ENDDO

                !--------------------------------------------------------------------------------------------
                ! Contribution to the induced dipole moment of the system.
                !--------------------------------------------------------------------------------------------
                dip_ind = dip_ind +  prefac * AIMAG(sk * MATMUL(sk_o, rk))
                !--------------------------------------------------------------------------------------------
                !--------------------------------------------------------------------------------------------

                !--------------------------------------------------------------------------------------------
                ! Contribution to the derivative of the induced dipole moment of the system.
                !--------------------------------------------------------------------------------------------
                IF (compute_der) THEN
                   rk_out = prefac * outer(rk, rk)
                   DO i = 1, nmol
                      iatom = 3 * (i-1) + 1
                      !Sum for i-th oxygen atom.
                      dip_ind_der(iatom, :, :) = dip_ind_der(iatom, :, :) + MATMUL(alpha(i, :, :), rk_out) * REAL(exp_ikr(i) * sk, KIND=KIND(1.0d0))
                      DO l = 1, 3
                         latom = iatom + l - 1
                         !Derivatives with respect to atom l of molecule j.
                         tnsr_tmp = REAL(sk_i(latom) * MATMUL(sk_o, rk_out), KIND=KIND(1.0d0))
                         CALL dip_ind_der_ij(dip_ind_der(iatom : iatom + 2, :, :), tnsr_tmp, l)
                         !Contribution from the derivative of alpha.
                         DO m = 1, 3
                            dip_ind_der(latom, m, :) = dip_ind_der(latom, m, :) + prefac * MATMUL(dalpha_dr(latom, :, m, :), rk) * AIMAG(exp_ikr(i) * sk)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
                !--------------------------------------------------------------------------------------------
                !--------------------------------------------------------------------------------------------

                !--------------------------------------------------------------------------------------------
                ! Contribution to the induced polarizability of the system.
                !--------------------------------------------------------------------------------------------
                sk_o_times_rk = MATMUL(sk_o, rk)
                pol_ind = pol_ind + prefac * REAL(outer_cmplx(sk_o_times_rk, CONJG(sk_o_times_rk)), KIND=KIND(1.0d0))
                !--------------------------------------------------------------------------------------------
                !--------------------------------------------------------------------------------------------
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !--------------------------------------------------------------------------------------------------
    ! Non-zero constant contribution to the induced polarizability of the system, from kx=ky=kz=0
    !--------------------------------------------------------------------------------------------------
    sum_alpha = SUM(alpha, DIM=1)
    pol_ind = pol_ind + f * MATMUL(sum_alpha, sum_alpha) / 3.0d0
    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------

  END SUBROUTINE long_range_ew

  SUBROUTINE dip_ind_der_ij(dmui_drj, tnsr, k)

    DOUBLE PRECISION, INTENT(INOUT) :: dmui_drj(3, 3, 3)
    DOUBLE PRECISION, INTENT(IN) :: tnsr(3, 3)
    INTEGER, INTENT(IN) :: k

    INTEGER :: h

    IF (k .EQ. 1) THEN
       !j-th molecule oxygen atom. Includes i==j term.
       dmui_drj(k, :, :) = dmui_drj(k, :, :) - gam * tnsr(:, :)
       !j-th molecule hydrogen atoms (part that comes from M site).
       DO h = 1, 2
          dmui_drj(k + h, :, :) = dmui_drj(k + h, :, :) - gam2 * tnsr(:, :)
       ENDDO
    ELSE
       !j-th molecule hydrogen atoms (part that comes from hydrogen atom positions).
       dmui_drj(k, :, :) = dmui_drj(k, :, :) - tnsr(:, :)
    ENDIF   

  END SUBROUTINE dip_ind_der_ij

  FUNCTION calc_dipole(atoms, charges, dip_ind)

    DOUBLE PRECISION, INTENT(IN) :: atoms(nat,3), charges(3), dip_ind(3)

    DOUBLE PRECISION :: calc_dipole(3)

    INTEGER :: i, iatom, k

    calc_dipole = 0.0d0
    !Permanent dipoles.
    DO i = 1, nmol
       iatom = 3 * (i - 1)
       DO k = 1, 3
          calc_dipole(:) = calc_dipole(:) + charges(k) * atoms(iatom + k, :) 
       ENDDO
       !Induced dipoles from electrostatic interaction with charges on other water molecules.
    ENDDO
    calc_dipole = calc_dipole + dip_ind

  END FUNCTION calc_dipole

  FUNCTION calc_dipole_derivative(charges, dip_ind_der) RESULT(dip_der)

    DOUBLE PRECISION, INTENT(IN) :: charges(3)
    DOUBLE PRECISION, INTENT(IN) :: dip_ind_der(nat, 3, 3)

    DOUBLE PRECISION :: dip_der(nat, 3, 3)

    INTEGER :: iatom, i, k

    dip_der = 0.0d0
    DO i=1, nmol
       iatom = 3 * (i - 1) + 1
       DO k = 1, 3
          !Gradient of the permanent molecular dipole moment.
          dip_der(iatom, k, k) = dip_der(iatom, k, k) + gam * charges(1) !i-th oxygen.
          dip_der(iatom+1:iatom+2, k, k) = dip_der(iatom+1:iatom+2, k, k) + charges(2) + gam2 * charges(1) !i-th hydrogen.
       ENDDO
    ENDDO
    dip_der = dip_der + dip_ind_der

  END FUNCTION calc_dipole_derivative

  FUNCTION calc_polarizability(pol_ind) RESULT(pol)

    DOUBLE PRECISION, INTENT(IN) :: pol_ind(3,3)

    DOUBLE PRECISION :: pol(3, 3)

    pol = 0.0d0
    pol = SUM(alpha, DIM=1)
    pol = pol - pol_ind 

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

    DOUBLE PRECISION, INTENT(IN) :: a(3), b(3)

    DOUBLE PRECISION :: outer(3,3)

    INTEGER :: i

    DO i=1, 3
       outer(i, :)  = a(i) * b(:)
    ENDDO

  END FUNCTION outer

  FUNCTION outer_cmplx(a, b)

    DOUBLE COMPLEX, INTENT(IN) :: a(3), b(3)

    DOUBLE COMPLEX :: outer_cmplx(3,3)

    INTEGER :: i

    DO i=1, 3
       outer_cmplx(i, :)  = a(i) * b(:)
    ENDDO

  END FUNCTION outer_cmplx

  FUNCTION cross_product(a, b) RESULT(c)

    DOUBLE PRECISION, INTENT(IN) :: a(3), b(3)

    DOUBLE PRECISION :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

  END FUNCTION cross_product

END SUBROUTINE h2o_dipole

