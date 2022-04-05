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
SUBROUTINE h2o_dipole(nat, atoms, dip, dip_der)

  IMPLICIT NONE
  
  INTEGER :: nat
  DOUBLE PRECISION :: atoms(nat, 3)
  DOUBLE PRECISION :: dip(3)
  DOUBLE PRECISION :: dip_der(nat, 3)
  
  DOUBLE PRECISION, PARAMETER :: qm = -1.1128d0
  DOUBLE PRECISION, PARAMETER :: qh = - qm / 2
  DOUBLE PRECISION, PARAMETER :: charges(3) = (/ qm, qh, qh /) !(n_charged_atom_per_mol)
  DOUBLE PRECISION, PARAMETER :: gam = 0.73612d0, gam2 = 0.5d0 * (1-gam)
  DOUBLE PRECISION, PARAMETER :: a_iso = 1.47d0 * 1.88973d0**3
  DOUBLE PRECISION :: alpha(3) !(xyz)
  
  DOUBLE PRECISION :: dip_i(nat/3, 3) !Molecular dipole moments. (nmol, xyz)
  DOUBLE PRECISION :: dip_der_i(nat/3, nat, 3, 3) !Gradients of molecular dipole moments. (nmol, natoms, xyz, xyz)

  DOUBLE PRECISION :: atoms_charged(nat, 3) !List of atoms, with O replaced by M site ('charged atoms'). (natoms, xyz)
  DOUBLE PRECISION :: E_stat(3) !(xyz)
  DOUBLE PRECISION :: r_ij(3), r, r2, r3
  DOUBLE PRECISION :: T_tnsr(nat/3, 3, 3, 3) !(nmol, n_charged_atom_per_mol, xyz, xyz)
  
  INTEGER :: i, iatom, j, jatom, k
  
  !Molecular polarizability.
  alpha = a_iso
  
  !Charged atoms (Site M instead of oxygen).
  atoms_charged = atoms
  DO i=1, nat, 3
      atoms_charged(i, :) = gam * atoms(i, :) + gam2 * (atoms(i+1, :) + atoms(i+2, :))
  ENDDO
  
  DO i=1, nat/3
     iatom = 3*(i-1) + 1
     !Coordinate-dependent quantities: E_stat and T_tnsr.
     CALL calc_E_stat_and_T(i, E_stat, T_tnsr)
     !Calculate dipole of i-th molecule. Using scaled charges for permanent dipole according to Hamm.
     dip_i(i, :) = calc_dipole(atoms_charged(iatom:iatom+2, :), charges / 1.3d0, E_stat)
     !Gradient of the i-th dipole moment.
     dip_der_i(i, :, : ,:) = calc_dipole_derivative(charges, T_tnsr, i)
  ENDDO
  
  dip = SUM(dip_i, DIM=1)
  dip_der = SUM(dip_der_i(:, :, 3, :), DIM=1)

 CONTAINS

  SUBROUTINE calc_E_stat_and_T(i, E_stat, T_tnsr)

    INTEGER :: i !Index of molecule for which we want to compute static electric field from surrounding molecules.
    DOUBLE PRECISION :: E_stat(3), T_tnsr(nat/3, 3, 3, 3)

    E_stat = 0.0d0
    T_tnsr = 0.0d0
    DO j = 1, nat/3
       IF (i == j) CYCLE
       DO k = 1, 3
          jatom = 3 * (j - 1) + k
          r_ij(:) = atoms(iatom, :) - atoms_charged(jatom, :)
          r2 = SUM(r_ij**2)
          r = SQRT(r2)
          r3 = r*r2
          E_stat = E_stat + charges(k) * r_ij / r3
          T_tnsr(j, k, :, :) = ( r2 - 3 * outer(r_ij, r_ij) ) / (r3 * r2)
       ENDDO
    ENDDO

  END SUBROUTINE calc_E_stat_and_T

  FUNCTION calc_dipole(atoms, charges, E_stat)

    DOUBLE PRECISION :: atoms(3,3), charges(3), E_stat(3)

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

  FUNCTION calc_dipole_derivative(charges, T_tnsr, i) RESULT(dip_der)

    DOUBLE PRECISION :: charges(3), T_tnsr(nat/3, 3, 3, 3)
    INTEGER :: i

    DOUBLE PRECISION :: dip_der(nat, 3, 3)

    INTEGER :: iatom, j, jatom, k

    dip_der = 0.0d0
    iatom = 3 * (i - 1) + 1
    DO j = 1, nat/3
       IF(i==j) CYCLE
       jatom = 3*(j-1) + 1
       !Sum for i-th oxygen atom.
       DO k = 1, 3
          dip_der(iatom, :, :) =  dip_der(iatom, :, :) + charges(k) * T_tnsr(j, k, :, :)
       ENDDO
       !j-th molecule oxygen atom.
       dip_der(jatom, :, :) = - gam * qm * T_tnsr(j, 1, :, :)
       !j-th molecule hydrogen atoms.
       DO k = 2, 3
          dip_der(jatom + k - 1, :, :) = - gam2 * qm * T_tnsr(j, 1, :, :) - qh * T_tnsr(j, k, :, :)
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

  FUNCTION outer(a, b)

    DOUBLE PRECISION :: a(3), b(3)

    DOUBLE PRECISION :: outer(3,3)

    INTEGER :: i

    DO i=1, 3
       outer(i, :)  = a(i) * b(:)
    ENDDO

  END FUNCTION outer

END SUBROUTINE h2o_dipole

