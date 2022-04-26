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
SUBROUTINE h2o_dipole(box, nat, atoms, compute_der, dip, dip_der)

  IMPLICIT NONE
  
  DOUBLE PRECISION, INTENT(IN) :: box(3)
  INTEGER, INTENT(IN) :: nat
  DOUBLE PRECISION, INTENT(IN) :: atoms(nat, 3)
  LOGICAL, INTENT(IN) :: compute_der
  DOUBLE PRECISION, INTENT(INOUT) :: dip(3)
  DOUBLE PRECISION, INTENT(INOUT) :: dip_der(nat, 3)

  DOUBLE PRECISION, PARAMETER :: pi = DACOS(-1.0d0), twopi = 2 * pi, sqrtpi = SQRT(pi)
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
  DOUBLE PRECISION :: dEdr(nat/3, nat, 3, 3) !(nmol, nat, xyz, xyz)
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
  
  !Coordinate-dependent quantities: E_stat and dEdr. This is where we worry about Ewald summation.
  CALL calc_el_field(E_stat, dEdr)

  DO i=1, nat/3
     iatom = 3*(i-1) + 1
     !Calculate dipole of i-th molecule. Using scaled charges for permanent dipole according to Hamm.
     dip_i(i, :) = calc_dipole(atoms_charged(iatom:iatom+2, :), charges / 1.3d0, E_stat(i, :))
     !Gradient of the i-th dipole moment.
     IF (compute_der) THEN
        dip_der_i(i, :, : ,:) = calc_dipole_derivative(charges / 1.3d0, dEdr(i, :, :, :), i)
     ELSE
        dip_der_i(i, :, :, :) = 0.0d0
     ENDIF
  ENDDO
  
  dip = SUM(dip_i, DIM=1)
  dip_der = SUM(dip_der_i(:, :, 3, :), DIM=1)

 CONTAINS

  SUBROUTINE calc_el_field(E_stat, dEdr)

    DOUBLE PRECISION, INTENT(INOUT) :: E_stat(nat/3, 3), dEdr(nat/3, nat, 3, 3)

    DOUBLE PRECISION :: a, rcut

    !----------------------------------------------
    ! Parameters for Ewald (calculated only once).
    !----------------------------------------------
    rcut = 1.5d0 * MINVAL(box) * MIN(0.5d0,1.2d0*nat**(-1.d0/6.d0))
    a = pi/rcut

    !----------------------------------------------
    ! Short-range part - sum over pairs.
    !----------------------------------------------
    CALL short_range_ew(atoms_charged, ro, a, rcut, E_stat, dEdr)

    !----------------------------------------------
    ! Long-range part - performs sum in k space.
    !----------------------------------------------
    CALL long_range_ew(atoms_charged, ro, a, E_stat, dEdr)

  END SUBROUTINE calc_el_field

  SUBROUTINE short_range_ew(r, ro, a, rcut, E_stat, dEdr)

    DOUBLE PRECISION, INTENT(IN) :: r(nat, 3), ro(nat/3, 3), a, rcut 
    DOUBLE PRECISION, INTENT(INOUT) :: E_stat(nat/3, 3), dEdr(nat/3, nat, 3, 3)

    INTEGER :: i, j, k, jatom, nx, ny, nz, nmax
    DOUBLE PRECISION :: r_ij_0(3), r_ij(3), dr, dr2, dr3, rcut2, screen, a_dr, gauss_part, T_tnsr(3, 3)
    LOGICAL :: self_term

    nmax=NINT(rcut/MINVAL(box))
    E_stat = 0.0d0
    dEdr = 0.0d0
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
                            T_tnsr = T_tnsr + charges(k) * short_range_T_tnsr(r_ij, dr2, dr3, a_dr, gauss_part, screen)
                         ENDIF
                         !--------------------------------------------------------------------------------------
                         ! Contribution to the electric field of molecule i from atom k in molecule j in cell n.
                         !--------------------------------------------------------------------------------------
                         E_stat(i, :) = E_stat(i, :) + charges(k) * r_ij(:) / dr3 * screen
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
                !Sum for i-th oxygen atom.
                dEdr(i, iatom, :, :) =  dEdr(i, iatom, :, :) + T_tnsr(:, :)
                !Derivative of electric field at molecule i w.r.t. coordinates of molecule j (includes self-term i==j).
                CALL el_field_der_ij(dEdr(i, jatom + 1: jatom + 3, :, :), T_tnsr(:, :), k)
             ENDIF
             !--------------------------------------------------------------------------------------------------------
             !--------------------------------------------------------------------------------------------------------
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

  SUBROUTINE long_range_ew(r, ro, a, E_stat, dEdr)

    DOUBLE PRECISION, INTENT(IN) :: r(nat, 3), ro(nat/3, 3), a    
    DOUBLE PRECISION, INTENT(INOUT) :: E_stat(nat/3, 3), dEdr(nat/3, nat, 3, 3)

    INTEGER :: i, j, k, jatom, kx, ky, kz, kmax
    DOUBLE PRECISION :: b, f, rk(3), rk_out(3, 3), rk2, rkmax2, lat(3), q(nat), re_part
    DOUBLE COMPLEX :: sk_i(nat), sk, exp_ikr(nat/3), tmp(nat/3)

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
             IF(compute_der) rk_out = outer(rk, rk)
             IF (rk2 .LT. rkmax2 .AND. rk2 .GT. EPSILON(0.d0)) THEN
                sk_i = f * (EXP(-b*rk2) / rk2) * q * EXP(-IU * k_dot_r(nat, rk, r))
                sk = SUM(sk_i)
                exp_ikr = EXP(IU * k_dot_r(nat/3, rk, ro))
                tmp(:) = exp_ikr(:) * SUM(sk_i)
                DO i = 1, nat/3
                   iatom = 3 * (i-1) + 1
                   E_stat(i,:) = E_stat(i, :) + rk(:)  * AIMAG(tmp(i))
                   IF (compute_der) THEN
                      !Sum for i-th oxygen atom.
                      dEdr(i, iatom, :, :) = dEdr(i, iatom, :, :) + rk_out(:, :) * REAL(tmp(i), KIND=KIND(1.0d0))
                      DO j = 1, nat/3
                         jatom = 3 * (j-1)
                         DO k = 1, 3
                            !Derivatives with respect to atom k of molecule j.
                            re_part = REAL(exp_ikr(i) * sk_i(jatom + k), KIND=KIND(1.0d0))
                            CALL el_field_der_ij(dEdr(i, jatom + 1: jatom + 3, :, :), re_part * rk_out(:, :), k)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE long_range_ew

  SUBROUTINE el_field_der_ij(dEi_drj, tnsr, k)

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

  END SUBROUTINE el_field_der_ij

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

  FUNCTION calc_dipole_derivative(charges, dEdr, i) RESULT(dip_der)

    DOUBLE PRECISION, INTENT(IN) :: charges(3)
    DOUBLE PRECISION, INTENT(IN) :: dEdr(nat, 3, 3)
    INTEGER, INTENT(IN) :: i

    DOUBLE PRECISION :: dip_der(nat, 3, 3)

    INTEGER :: iatom, j, k

    dip_der = 0.0d0
    iatom = 3 * (i - 1) + 1
    DO k = 1, 3
       !Gradient of the permanent molecular dipole moment.
       dip_der(iatom, k, k) = gam * charges(1) !i-th oxygen.
       dip_der(iatom+1:iatom+2, k, k) = charges(2) + gam2 * charges(1) !i-th hydrogen.
       !Multiply induced part of dipole derivative by alpha.
       DO j = 1, nat
          dip_der(j, k, :) = dip_der(j, k, :) + alpha(k) * dEdr(j, k, :)
       ENDDO
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

