! The main program which runs our driver test case potentials
!
! Copyright (C) 2013, Joshua More and Michele Ceriotti
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!
! Currently the potentials implemented are the Lennard-Jones
! potential, the Silvera-Goldman para-hydrogen potential and
! the ideal gas (i.e. no interaction at all)

      PROGRAM DRIVER
         USE LJ
         USE LJPolymer
         USE SG
         USE PSWATER
         USE F90SOCKETS, ONLY : open_socket, writebuffer, readbuffer
      IMPLICIT NONE
      
      ! SOCKET VARIABLES
      INTEGER, PARAMETER :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
      INTEGER socket, inet, port        ! socket ID & address of the server
      CHARACTER(LEN=1024) :: host
      
      ! COMMAND LINE PARSING
      CHARACTER(LEN=1024) :: cmdbuffer
      INTEGER ccmd, vstyle, vseed
      INTEGER, ALLOCATABLE :: seed(:)
      INTEGER verbose
      INTEGER commas(4), par_count      ! stores the index of commas in the parameter string
      DOUBLE PRECISION vpars(5)         ! array to store the parameters of the potential
      
      ! SOCKET COMMUNICATION BUFFERS
      CHARACTER(LEN=12) :: header
      LOGICAL :: isinit=.false., hasdata=.false.
      INTEGER cbuf, rid, length
      CHARACTER(LEN=4096) :: initbuffer      ! it's unlikely a string this large will ever be passed...
      CHARACTER(LEN=4096) :: string,string2,trimmed  ! it's unlikely a string this large will ever be passed...
      DOUBLE PRECISION, ALLOCATABLE :: msgbuffer(:)
      
      ! PARAMETERS OF THE SYSTEM (CELL, ATOM POSITIONS, ...)
      DOUBLE PRECISION sigma, eps, rc, rn, ks ! potential parameters
      DOUBLE PRECISION stiffness ! lennard-jones polymer
      INTEGER n_monomer ! lennard-jones polymer
      INTEGER nat
      DOUBLE PRECISION pot, dpot, dist
      DOUBLE PRECISION, ALLOCATABLE :: atoms(:,:), forces(:,:), datoms(:,:)
      DOUBLE PRECISION cell_h(3,3), cell_ih(3,3), virial(3,3), mtxbuf(9), dip(3), charges(3), dummy(3,3,3), vecdiff(3)
      DOUBLE PRECISION, ALLOCATABLE :: friction(:,:)
      DOUBLE PRECISION volume
      DOUBLE PRECISION, PARAMETER :: fddx = 1.0d-5
      
      ! NEIGHBOUR LIST ARRAYS
      INTEGER, DIMENSION(:), ALLOCATABLE :: n_list, index_list
      DOUBLE PRECISION init_volume, init_rc ! needed to correctly adjust the cut-off radius for variable cell dynamics
      DOUBLE PRECISION, ALLOCATABLE :: last_atoms(:,:) ! Holds the positions when the neighbour list is created
      DOUBLE PRECISION displacement ! Tracks how far each atom has moved since the last call of nearest_neighbours

      ! DMW
      DOUBLE PRECISION efield(3)
      INTEGER i, j

      !------------------------------------------------------------------
      !------------------------------------------------------------------
      ! Inserted code by TB:
      DOUBLE PRECISION :: rnorm, r(3), r2, r5, r7, pol ! Interatomic distance, difference between their z components and polarizability.
      DOUBLE PRECISION, ALLOCATABLE :: pol_der(:, :)
      DOUBLE PRECISION, ALLOCATABLE :: dip_der(:, :)
      DOUBLE PRECISION :: polar(3, 3)
      DOUBLE PRECISION :: cell_abc(3)
      DOUBLE PRECISION, PARAMETER :: atomic_pol = 2.67234d0 ! Atomic polarizability in bohr^3.
      CHARACTER(LEN=1000000) :: out_string ! it's unlikely a string this large will ever be passed...
      INTEGER :: nbeads=1, dipole_freq=1, der_freq=1, nsteps_eq=1000000, nsteps_neq=1000000
      INTEGER :: counter, step
      LOGICAL :: compute_dip, compute_der, neq = .FALSE.
      ! Used for numerical derivative of the dipole moment, only for testing.
      !DOUBLE PRECISION, ALLOCATABLE :: dip_der_num(:, :)
      !DOUBLE PRECISION :: dip_tmp(3)
      !DOUBLE PRECISION, ALLOCATABLE :: atoms_tmp(:, :)

 
      !------------------------------------------------------------------
      !------------------------------------------------------------------
      
      ! parse the command line parameters
      ! intialize defaults
      ccmd = 0
      inet = 1
      host = "localhost"//achar(0)
      port = 31415
      verbose = 0
      par_count = 0
      vstyle = -1
      rc = 0.0d0
      init_rc = 0.0d0
      volume = 0.0d0
      init_volume = 0.0d0

      DO i = 1, COMMAND_ARGUMENT_COUNT()
         CALL GET_COMMAND_ARGUMENT(i, cmdbuffer)
         IF (cmdbuffer == "-u") THEN ! flag for unix socket
            inet = 0
            ccmd = 0
         ELSEIF (cmdbuffer == "-h") THEN ! read the hostname (deprecated)
            ccmd = 1
         ELSEIF (cmdbuffer == "-a") THEN ! read the hostname (address)
            ccmd = 1
         ELSEIF (cmdbuffer == "-p") THEN ! reads the port number
            ccmd = 2
         ELSEIF (cmdbuffer == "-m") THEN ! reads the style of the potential function
            ccmd = 3
         ELSEIF (cmdbuffer == "-o") THEN ! reads the parameters
            ccmd = 4
         ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
            verbose = 1
         ELSEIF (cmdbuffer == "-vv") THEN ! flag for verbose standard output
            verbose = 2
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) " Unrecognized command line argument", ccmd
               CALL helpmessage
               STOP "ENDED"
            ENDIF
            IF (ccmd == 1) THEN
               host = trim(cmdbuffer)//achar(0)
            ELSEIF (ccmd == 2) THEN
               READ(cmdbuffer,*) port
            ELSEIF (ccmd == 3) THEN
               IF (trim(cmdbuffer) == "lj") THEN
                  vstyle = 1
               ELSEIF (trim(cmdbuffer) == "sg") THEN
                  vstyle = 2
               ELSEIF (trim(cmdbuffer) == "harm") THEN
                  vstyle = 3
               ELSEIF (trim(cmdbuffer) == "harm3d") THEN
                  vstyle = 30
               ELSEIF (trim(cmdbuffer) == "morse") THEN
                  vstyle = 4
               ELSEIF (trim(cmdbuffer) == "zundel") THEN
                  vstyle = 5
               ELSEIF (trim(cmdbuffer) == "qtip4pf") THEN
                  vstyle = 6
               ELSEIF (trim(cmdbuffer) == "linear") THEN
                  vstyle = 7
               ELSEIF (trim(cmdbuffer) == "pswater") THEN
                  vstyle = 8
               ELSEIF (trim(cmdbuffer) == "lepsm1") THEN
                  vstyle = 9
               ELSEIF (trim(cmdbuffer) == "lepsm2") THEN
                  vstyle = 10
               ELSEIF (trim(cmdbuffer) == "qtip4pf-efield") THEN
                  vstyle = 11
               ELSEIF (trim(cmdbuffer) == "eckart") THEN
                  vstyle = 20
               ELSEIF (trim(cmdbuffer) == "ch4hcbe") THEN
                  vstyle = 21
               ELSEIF (trim(cmdbuffer) == "ljpolymer") THEN
                  vstyle = 22
               ELSEIF (trim(cmdbuffer) == "MB") THEN
                  vstyle = 23
               ELSEIF (trim(cmdbuffer) == "doublewell") THEN
                  vstyle = 25
               ELSEIF (trim(cmdbuffer) == "doublewell_1D") THEN
                  vstyle = 24
               ELSEIF (trim(cmdbuffer) == "water_dip_pol") THEN
                  vstyle = 26
               ELSEIF (trim(cmdbuffer) == "gas") THEN
                  vstyle = 0  ! ideal gas
               ELSEIF (trim(cmdbuffer) == "dummy") THEN
                  vstyle = 99 ! returns non-zero but otherwise meaningless values
               ELSE
                  WRITE(*,*) " Unrecognized potential type ", trim(cmdbuffer)
                  WRITE(*,*) " Use -m [dummy|gas|lj|sg|harm|harm3d|morse|zundel|qtip4pf|pswater|lepsm1|lepsm2|qtip4pf-efield|eckart|ch4hcbe|ljpolymer|MB|doublewell|doublewell_1D] "
                  STOP "ENDED"
               ENDIF
            ELSEIF (ccmd == 4) THEN
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vpars(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vpars(par_count)
            ENDIF
            ccmd = 0
         ENDIF
      ENDDO
      
      IF (vstyle == -1) THEN
         WRITE(*,*) " Error, type of potential not specified."
         CALL helpmessage
         STOP "ENDED"
      ELSEIF (0 == vstyle) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for ideal gas."
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (99 == vstyle) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for dummy output."
            STOP "ENDED"
         ENDIF
         CALL RANDOM_SEED(size=vseed)
         ALLOCATE(seed(vseed))
         seed = 12345
         CALL RANDOM_SEED(put=seed)
         isinit = .true.         
      ELSEIF (6 == vstyle) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error:  no initialization string needed for qtip4pf."
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (11== vstyle) THEN
         IF (par_count .ne. 3) THEN
            WRITE(*,*) "Error:  incorrect initialization string included for qtip4pf-efield. &
     &             Provide the three components of the electric field in V/nm"
            STOP "ENDED"
         ELSE
            ! We take in an electric field in volts / nm.This must be converted to Eh / (e a0).
            do i=1,3
             efield(i) = vpars(i) / 5.14220652d2
            enddo
         ENDIF
         isinit = .true.
      ELSEIF (5 == vstyle) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for zundel."
            STOP "ENDED"
         ENDIF
         CALL prezundelpot()
         CALL prezundeldip()
         isinit = .true.
      ELSEIF (21 == vstyle) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for CH4+H CBE potential."
            STOP "ENDED"
         ENDIF
         CALL prepot()
         isinit = .true.
      ELSEIF (4 == vstyle) THEN
         IF (par_count == 0) THEN ! defaults (OH stretch)
            vpars(1) = 1.8323926 ! r0
            vpars(2) = 0.18748511263179304 ! D
            vpars(3) = 1.1562696428501682 ! a
         ELSEIF ( 2/= par_count) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For morse potential use -o r0,D,a (in a.u.) "
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (20 == vstyle) THEN !eckart
         IF (par_count == 0) THEN ! defaults values 
            vpars(1) = 0.0d0
            vpars(2) = 0.66047 
            vpars(3) = (6*12)/( 1836 * (vpars(2)**2) *( (4.D0 * ATAN(1.0d0) )**2 ) )
            vpars(4) = 1836*(3800.0d0/219323d0)**2
         ELSEIF ( 4/= par_count) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For eckart potential use  AA,A,B,k"
            STOP "ENDED"
         ENDIF
         isinit = .true.

      ELSEIF (23 == vstyle) THEN !MB
         IF (par_count == 0) THEN ! defaults values 
            vpars(1) = 0.004737803248674678
         ELSEIF ( 1/= par_count) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For MB potential up to 1 param can be specified"
            STOP "ENDED"
         ENDIF
         isinit = .true.
      ELSEIF (22 == vstyle) THEN !ljpolymer
         IF (4/= par_count) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For ljpolymer potential use n_monomer,sigma,epsilon,cutoff"
            STOP "ENDED"
         ELSE
            n_monomer = nint(vpars(1))
            sigma = vpars(2)
            eps = vpars(3)
            rc = vpars(4)
            rn = rc * 1.2d0
            stiffness = 36.d0 * (2.d0 ** (2.d0/3.d0))*eps
            isinit = .true.
         ENDIF
      ELSEIF (vstyle == 8) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for Partridge-Schwenke H2O potential."
            STOP "ENDED"
         END IF
      ELSEIF (vstyle == 9) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for LEPSM1."
            STOP "ENDED"
         END IF
      ELSEIF (vstyle == 10) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for LEPSM2."
            STOP "ENDED" 
         ENDIF   
         isinit = .true.
      ELSEIF (vstyle == 11) THEN
         IF (par_count .ne. 3) THEN
            WRITE(*,*) "Error:  incorrect initialization string included for qtip4pf-efield. &
     &    Provide the three components of the electric field in V/nm"
            STOP "ENDED"
         ELSE
            ! We take in an electric field in volts / nm.This must be converted 
            ! to Eh / (e a0).
            do i=1,3
             efield(i) = vpars(i) / 5.14220652d2
            enddo
         ENDIF
         isinit = .true.
      ELSEIF (vstyle == 1) THEN
         IF (par_count /= 3) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For LJ potential use -o sigma,epsilon,cutoff "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         sigma = vpars(1)
         eps = vpars(2)
         rc = vpars(3)
         rn = rc*1.2
         isinit = .true.
      ELSEIF (vstyle == 2) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For SG potential use -o cutoff "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         rc = vpars(1)
         rn = rc*1.2
         isinit = .true.
      ELSEIF (vstyle == 3) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For 1D harmonic potential use -o k "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         ks = vpars(1)
         isinit = .true.
      ELSEIF (vstyle == 30) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For 3D harmonic potential use -o k "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         ks = vpars(1)  !è la k dell'ho, unica chiaramente perché in 1D
         isinit = .true.
      ELSEIF (vstyle == 7) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For a linear potential use -o k "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         ks = vpars(1)
         isinit = .true.

      ELSEIF (25 == vstyle) THEN !doublewell
         IF ( par_count /= 0 ) THEN
                 WRITE(*,*) "Error: no initialization string needed for doublewell."
            STOP "ENDED" 
         ENDIF   
         isinit = .true.

      ELSEIF (24 == vstyle) THEN !doublewell_1D
         IF ( par_count /= 0 ) THEN
                 WRITE(*,*) "Error: no initialization string needed for 1-dimensional doublewell."
            STOP "ENDED" 
         ENDIF   
         isinit = .true.
      ELSEIF (vstyle == 26) THEN !water dipole and polarizability
         counter = 0
         IF (par_count == 5) THEN
            nbeads=INT(vpars(1))
            dipole_freq=INT(vpars(2))
            der_freq=INT(vpars(3))
            nsteps_eq=INT(vpars(4))
            nsteps_neq=INT(vpars(5))
         ENDIF
         isinit = .true.
      ENDIF

      IF (verbose > 0) THEN
         WRITE(*,*) " DRIVER - Connecting to host ", trim(host)
         IF (inet > 0) THEN
            WRITE(*,*) " on port ", port, " using an internet socket."
         ELSE
            WRITE(*,*) " using an UNIX socket."
         ENDIF
      ENDIF

      ! Calls the interface to the POSIX sockets library to open a communication channel
      CALL open_socket(socket, inet, port, host)
      nat = -1
      DO WHILE (.true.) ! Loops forever (or until the wrapper ends!)

         ! Reads from the socket one message header
         CALL readbuffer(socket, header, MSGLEN)
         IF (verbose > 0) WRITE(*,*) " Message from server: ", trim(header)

         IF (trim(header) == "STATUS") THEN
            ! The wrapper is inquiring on what we are doing
            IF (.not. isinit) THEN
               CALL writebuffer(socket,"NEEDINIT    ",MSGLEN)  ! Signals that we need initialization data
               IF (verbose > 1) WRITE(*,*) "    !write!=> ", "NEEDINIT    "
            ELSEIF (hasdata) THEN
               CALL writebuffer(socket,"HAVEDATA    ",MSGLEN)  ! Signals that we are done computing and can return forces
               IF (verbose > 1) WRITE(*,*) "    !write!=> ", "HAVEDATA    "
            ELSE
               CALL writebuffer(socket,"READY       ",MSGLEN)  ! We are idling and eager to compute something
               IF (verbose > 1) WRITE(*,*) "    !write!=> ", "READY       "
            ENDIF
         ELSEIF (trim(header) == "INIT") THEN     ! The driver is kindly providing a string for initialization
            CALL readbuffer(socket, rid)
            IF (verbose > 1) WRITE(*,*) "    !read!=> RID: ", rid
            CALL readbuffer(socket, cbuf)
            IF (verbose > 1) WRITE(*,*) "    !read!=> init_length: ", cbuf
            CALL readbuffer(socket, initbuffer, cbuf)
            IF (verbose > 1) WRITE(*,*) "    !read!=> init_string: ", cbuf
            IF (verbose > 0) WRITE(*,*) " Initializing system from wrapper, using ", trim(initbuffer)
            isinit=.true. ! We actually do nothing with this string, thanks anyway. Could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver
         ELSEIF (trim(header) == "POSDATA") THEN  ! The driver is sending the positions of the atoms. Here is where we do the calculation!

            ! Parses the flow of data from the socket
            CALL readbuffer(socket, mtxbuf, 9)  ! Cell matrix
            IF (verbose > 1) WRITE(*,*) "    !read!=> cell: ", mtxbuf
            cell_h = RESHAPE(mtxbuf, (/3,3/))
            CALL readbuffer(socket, mtxbuf, 9)  ! Inverse of the cell matrix (so we don't have to invert it every time here)
            IF (verbose > 1) WRITE(*,*) "    !read!=> cell-1: ", mtxbuf
            cell_ih = RESHAPE(mtxbuf, (/3,3/))

            ! The wrapper uses atomic units for everything, and row major storage.
            ! At this stage one should take care that everything is converted in the
            ! units and storage mode used in the driver.
            cell_h = transpose(cell_h)
            cell_ih = transpose(cell_ih)
            ! We assume an upper triangular cell-vector matrix
            volume = cell_h(1,1)*cell_h(2,2)*cell_h(3,3)

            CALL readbuffer(socket, cbuf)       ! The number of atoms in the cell
            IF (verbose > 1) WRITE(*,*) "    !read!=> cbuf: ", cbuf
            IF (nat < 0) THEN  ! Assumes that the number of atoms does not change throughout a simulation, so only does this once
               nat = cbuf
               IF (verbose > 0) WRITE(*,*) " Allocating buffer and data arrays, with ", nat, " atoms"
               ALLOCATE(msgbuffer(3*nat))
               ALLOCATE(atoms(nat,3), datoms(nat,3))
               ALLOCATE(forces(nat,3))
               ALLOCATE(friction(3*nat,3*nat))
               atoms = 0.0d0
               datoms = 0.0d0
               forces = 0.0d0
               friction = 0.0d0
               msgbuffer = 0.0d0
            ENDIF

            CALL readbuffer(socket, msgbuffer, nat*3)
            IF (verbose > 1) WRITE(*,*) "    !read!=> positions: ", msgbuffer
            DO i = 1, nat
               atoms(i,:) = msgbuffer(3*(i-1)+1:3*i)
            ENDDO

            IF (vstyle == 0) THEN   ! ideal gas, so no calculation done
               pot = 0
               forces = 0.0d0
               virial = 1.0d-200   
               ! returns a tiny but non-zero stress, so it can
               ! bypass the check for zero virial that is used
               ! to avoid running constant-pressure simulations
               ! with a code that cannot compute the virial
            ELSEIF (vstyle == 99) THEN ! dummy output, useful to test that i-PI "just runs"
               call random_number(pot)
               pot = pot - 0.5                
               call random_number(forces)
               forces = forces - 0.5
               call random_number(virial)
               virial = virial - 0.5
            ELSEIF (vstyle == 3) THEN ! 1D harmonic potential, so only uses the first position variable
               pot = 0.5*ks*atoms(1,1)**2
               forces = 0.0d0
               forces(1,1) = -ks*atoms(1,1)
               virial = 0.0d0
               virial(1,1) = forces(1,1)*atoms(1,1)
               IF (.NOT. ALLOCATED(dip_der)) ALLOCATE (dip_der(nat, 3))
               dip = 0.0d0; dip(1) = atoms(1,1)
               dip_der(:,:) = 0.0d0; dip_der(1,1) = 1.0d0
               polar(:,:) = 0.0d0; polar(1,1) = atoms(1,1)**2 
            ELSEIF (vstyle == 30) THEN ! 3D harmonic potential
                   pot = 0.0d0
                   forces = 0.0d0
                   virial = 0.0d0
               DO i=1,nat 
                   pot = pot + 0.5*ks*atoms(i,1)**2 + 0.5*ks*atoms(i,2)**2 + 0.5*ks*atoms(i,3)**2
                   forces(i,1) = -ks*atoms(i,1)
                   forces(i,2) = -ks*atoms(i,2)
                   forces(i,3) = -ks*atoms(i,3)
                   virial(i,1) = forces(i,1)*atoms(i,1)
                   virial(i,2) = forces(i,2)*atoms(i,2)
                   virial(i,3) = forces(i,3)*atoms(i,3)
                enddo
            ELSEIF (vstyle == 7) THEN ! linear potential in x position of the 1st atom
               pot = ks*atoms(1,1)
               forces = 0.0d0
               virial = 0.0d0
               forces(1,1) = -ks
               virial(1,1) = forces(1,1)*atoms(1,1)
            ELSEIF (vstyle == 4) THEN ! Morse potential.
               IF (nat/=1) THEN
                  WRITE(*,*) "Expecting 1 atom for 3D Morse (use the effective mass for the atom mass to get proper frequency!) "
                  STOP "ENDED"
               ENDIF
               CALL getmorse(vpars(1), vpars(2), vpars(3), atoms, pot, forces)
            ELSEIF (vstyle == 5) THEN ! Zundel potential.
               IF (nat/=7) THEN
                  WRITE(*,*) "Expecting 7 atoms for Zundel potential, O O H H H H H "
                  STOP "ENDED"
               ENDIF

               CALL zundelpot(pot,atoms)
               CALL zundeldip(dip,atoms)

               datoms=atoms
               DO i=1,7  ! forces by finite differences
                  DO j=1,3
                     datoms(i,j)=atoms(i,j)+fddx
                     CALL zundelpot(dpot, datoms)
                     datoms(i,j)=atoms(i,j)-fddx
                     CALL zundelpot(forces(i,j), datoms)
                     datoms(i,j)=atoms(i,j)
                     forces(i,j)=(forces(i,j)-dpot)/(2*fddx)
                  ENDDO
               ENDDO
               ! do not compute the virial term

           ELSEIF (vstyle == 21) THEN ! CBE CH4+H potential.
               IF (nat/=6) THEN
                  WRITE(*,*) "Expecting 6 atoms for CH4+H potential, H, C, H, H, H, H "
                  WRITE(*,*) "The expected order is such that atoms 1 to 5 are reactant_1 (CH4)"
                  WRITE(*,*) "and atom 6 is reactant_2 ( H 'free') "
                  STOP "ENDED"
               ENDIF

               CALL ch4hpot_inter(atoms, pot)
               datoms=atoms
               DO i=1,6  ! forces by finite differences
                  DO j=1,3
                     datoms(i,j)=atoms(i,j)+fddx
                     CALL ch4hpot_inter(datoms, dpot)
                     datoms(i,j)=atoms(i,j)-fddx
                     CALL ch4hpot_inter(datoms, forces(i,j))
                     datoms(i,j)=atoms(i,j)
                     forces(i,j)=(forces(i,j)-dpot)/(2*fddx)
                  ENDDO
               ENDDO
               ! do not compute the virial term

            ELSEIF (vstyle == 6) THEN ! qtip4pf potential.
               IF (mod(nat,3)/=0) THEN
                  WRITE(*,*) " Expecting water molecules O H H O H H O H H but got ", nat, "atoms"
                  STOP "ENDED"
               ENDIF
               vpars(1) = cell_h(1,1)
               vpars(2) = cell_h(2,2)
               vpars(3) = cell_h(3,3)
               IF (cell_h(1,2).gt.1d-10 .or. cell_h(1,3).gt.1d-12  .or. cell_h(2,3).gt.1d-12) THEN                       
                  WRITE(*,*) " qtip4pf PES only works with orthorhombic cells", cell_h(1,2), cell_h(1,3), cell_h(2,3)
                  STOP "ENDED"
               ENDIF
               CALL qtip4pf(vpars(1:3),atoms,nat,forces,pot,virial)
               dip(:) = 0.0
               DO i=1, nat, 3
                  dip = dip -1.1128d0 * atoms(i,:) + 0.5564d0 * (atoms(i+1,:) + atoms(i+2,:))
               ENDDO
               ! do not compute the virial term
            ELSEIF (vstyle == 11) THEN ! efield potential.             
               IF (mod(nat,3)/=0) THEN
                  WRITE(*,*) " Expecting water molecules O H H O H H O H H but got ", nat, "atoms"
                  STOP "ENDED"
               ENDIF
               CALL efield_v(atoms,nat,forces,pot,virial,efield)
            ELSEIF (vstyle == 8) THEN ! PS water potential.
               IF (nat/=3) THEN
                  WRITE(*,*) "Expecting 3 atoms for P-S water potential, O H H "
                  STOP "ENDED"
               ENDIF

               dip=0.0
               vecdiff=0.0
               ! lets fold the atom positions back to center in case the water travelled far away. 
               ! this avoids problems if the water is splic across (virtual) periodic boundaries
               ! OH_1
               call vector_separation(cell_h, cell_ih, atoms(2,:), atoms(1,:), vecdiff, dist)
               atoms(2,:)=vecdiff(:)
               ! OH_2
               call vector_separation(cell_h, cell_ih, atoms(3,:), atoms(1,:), vecdiff, dist)
               atoms(3,:)=vecdiff(:)
               ! O in center
               atoms(1,:)=0.d0



               atoms = atoms*0.52917721d0    ! pot_nasa wants angstrom
               call pot_nasa(atoms, forces, pot)
               call dms_nasa(atoms, charges, dummy) ! MR: trying to print out the right charges
               dip(:)=atoms(1,:)*charges(1)+atoms(2,:)*charges(2)+atoms(3,:)*charges(3)
               ! MR: the above line looks like it provides correct results in eAngstrom for dipole! 
               pot = pot*0.0015946679     ! pot_nasa gives kcal/mol
               forces = forces * (-0.00084329756) ! pot_nasa gives V in kcal/mol/angstrom
               ! do not compute the virial term
            ELSEIF (vstyle == 9) THEN
               IF (nat /= 3) THEN
                  WRITE(*,*) "Expecting 3 atoms for LEPS Model 1  potential, A B C "
                  STOP "ENDED"
               END IF
               CALL LEPS_M1(3, atoms, pot, forces)
            ELSEIF (vstyle == 10) THEN
               IF (nat /= 3) THEN
                  WRITE(*,*) "Expecting 4 atoms for LEPS Model 2  potential, A B C D n"
                  STOP "ENDED"
               END IF
               CALL LEPS_M2(4, atoms, pot, forces)
               
            ELSEIF (vstyle == 20) THEN ! eckart potential.
               CALL geteckart(nat,vpars(1), vpars(2), vpars(3),vpars(4), atoms, pot, forces)

            ELSEIF (vstyle == 23) THEN ! MB.
               IF (nat/=1) THEN
                  WRITE(*,*) "Expecting 1 atom for MB"
                  STOP "ENDED"
               ENDIF
               !atoms = atoms*0.52917721d0  !Change to angstrom
               CALL get_MB(nat,vpars(1), atoms, pot, forces)
            ELSEIF (vstyle == 25) THEN ! qQ
               CALL getdoublewell(nat, atoms, pot, forces)
               CALL dw_friction(nat, atoms, friction)

            ELSEIF (vstyle == 24) THEN ! qQ
               CALL getdoublewell_1D(nat, atoms, pot, forces)
               CALL dw1d_friction(nat, atoms, friction)
               CALL dw1d_dipole(nat, atoms, dip)

            ELSEIF (vstyle == 26) THEN   ! Sets force and potential to zero,
                                         ! computes only dipole moment, its gradient, and polarizability.
               pot = 0
               forces = 0.0d0
               virial = 0.0d0 
               DO i = 1, 3
                  vpars(i) = cell_h(i, i)
               ENDDO
               IF (.NOT. ALLOCATED(dip_der)) ALLOCATE (dip_der(nat, 3))
               step = counter/nbeads
               IF ((step .GT. nsteps_eq) .AND. (.NOT. neq) ) THEN
                  !Step greater than number of equilibrium dynamics steps, so we mark that we
                  !are now doing neq dynamics and reset counter/step to zero.
                  neq = .TRUE.
                  counter = 0; step = 0
               ELSEIF ((step .GT. nsteps_neq) .AND. neq) THEN
                  !Step greater than number of nonequilibrium dynamics steps, so we reset counter/step to zero.
                  counter = 0; step = 0
               ENDIF
               compute_dip = MOD(step, dipole_freq) .EQ. 0                  !Check if we should compute dipole (and polarizability).
               compute_der = (MOD(step, der_freq) .EQ. 0) .AND. (.NOT. neq) !Check if we should compute dipole derivative.
               IF (compute_dip .OR. compute_der) THEN
                  CALL h2o_dipole(vpars(1:3), nat, atoms, compute_der, dip, dip_der, polar)
               ELSE
                  dip = 0; dip_der = 0; polar = 0
               ENDIF
               !Numerical derivative (computed at first step only): Gives same result as the analytical derivative.
               !Used for testing only.
               !IF ((step .EQ. 0) .AND. (.NOT. neq)) THEN 
               !IF (.NOT. ALLOCATED(dip_der_num)) ALLOCATE(dip_der_num(nat, 3))
               !   DO i = 1, nat
               !      DO j = 1, 3
               !         atoms_tmp = atoms
               !         atoms_tmp(i, j) = atoms(i, j) + 0.00005d0
               !         CALL h2o_dipole(vpars(1:3), nat, atoms_tmp, .FALSE., dip, dip_der, polar)
               !         dip_tmp = dip
               !         atoms_tmp(i, j) = atoms(i, j) - 0.00005d0
               !         CALL h2o_dipole(vpars(1:3), nat, atoms_tmp, .FALSE., dip, dip_der, polar)
               !         dip_der_num(i, j) = (dip_tmp(3) - dip(3)) / 0.0001d0
               !      ENDDO
               !   ENDDO
               !   dip_der = dip_der_num
               !ENDIF

               counter = counter + 1

            ELSE
               IF ((allocated(n_list) .neqv. .true.)) THEN
                  IF (verbose > 0) WRITE(*,*) " Allocating neighbour lists."
                  ALLOCATE(n_list(nat*(nat-1)/2))
                  ALLOCATE(index_list(nat))
                  ALLOCATE(last_atoms(nat,3))
                  last_atoms = 0.0d0
                  CALL nearest_neighbours(rn, nat, atoms, cell_h, cell_ih, index_list, n_list)
                  last_atoms = atoms
                  init_volume = volume
                  init_rc = rc
               ENDIF

               ! Checking to see if we need to re-calculate the neighbour list
               rc = init_rc*(volume/init_volume)**(1.0/3.0)
               DO i = 1, nat
                  CALL separation(cell_h, cell_ih, atoms(i,:), last_atoms(i,:), displacement)
                  ! Note that displacement is the square of the distance moved by atom i since the last time the neighbour list was created.
                  IF (4*displacement > (rn-rc)*(rn-rc)) THEN
                     IF (verbose > 0) WRITE(*,*) " Recalculating neighbour lists"
                     CALL nearest_neighbours(rn, nat, atoms, cell_h, cell_ih, index_list, n_list)
                     last_atoms = atoms
                     rn = 1.2*rc
                     EXIT
                  ENDIF
               ENDDO

               IF (vstyle == 1) THEN
                  CALL LJ_getall(rc, sigma, eps, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
                  DO i = 1, 3
                     cell_abc(i) = cell_h(i,i)
                  ENDDO
                  pol = 0.0d0
                  IF (.NOT. ALLOCATED(pol_der)) ALLOCATE (pol_der(nat, 3))
                  pol_der = 0.0d0
                  DO i=1, nat
                     DO j=1, nat
                        IF (j .ne. i) THEN
                           r = atoms(i, :) - atoms(j, :)
                           WHERE(r .GT. cell_abc/2) r = r - cell_abc
                           WHERE(r .LE. -cell_abc/2) r = r + cell_abc
                           r2 = SUM(r**2)
                           rnorm = SQRT(r2)
                           !IF (rnorm < 4*rc) THEN
                              r5 = rnorm**5
                              r7 = r5 * r2
                              pol = pol + (3 * r(3)**2 - r2) / r5
                              pol_der(i,1:2) = pol_der(i,1:2) + r(1:2) * (r2 - 5 * r(3) ** 2) / r7
                              pol_der(i,3) = pol_der(i,3) + r(3) * (3 * r2 - 5 * r(3) ** 2) / r7
                           !END IF
                        END IF
                     ENDDO
                  ENDDO
                  !pol = atomic_pol * (nat + pol * atomic_pol)
                  pol = atomic_pol * pol * atomic_pol
                  pol_der = 6 * atomic_pol**2 * pol_der
               ELSEIF (vstyle == 2) THEN
                  CALL SG_getall(rc, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
               ELSEIF (vstyle == 22) THEN ! ljpolymer potential.
                  CALL ljpolymer_getall(n_monomer,rc,sigma,eps,stiffness,nat,atoms,cell_h,cell_ih,index_list,n_list,pot,forces,virial)
               ENDIF
               IF (verbose > 0) WRITE(*,*) " Calculated energy is ", pot
            ENDIF
            hasdata = .true. ! Signal that we have data ready to be passed back to the wrapper
         ELSEIF (trim(header) == "GETFORCE") THEN  ! The driver calculation is finished, it's time to send the results back to the wrapper

            ! Data must be re-formatted (and units converted) in the units and shapes used in the wrapper
            DO i = 1, nat
               msgbuffer(3*(i-1)+1:3*i) = forces(i,:)
            ENDDO
            virial = transpose(virial)

            CALL writebuffer(socket,"FORCEREADY  ",MSGLEN)
            IF (verbose > 1) WRITE(*,*) "    !write!=> ", "FORCEREADY  "
            CALL writebuffer(socket,pot)  ! Writing the potential
            IF (verbose > 1) WRITE(*,*) "    !write!=> pot: ", pot
            CALL writebuffer(socket,nat)  ! Writing the number of atoms
            IF (verbose > 1) WRITE(*,*) "    !write!=> nat:", nat
            CALL writebuffer(socket,msgbuffer,3*nat) ! Writing the forces
            IF (verbose > 1) WRITE(*,*) "    !write!=> forces:", msgbuffer
            CALL writebuffer(socket,reshape(virial,(/9/)),9)  ! Writing the virial tensor, NOT divided by the volume
            IF (verbose > 1) WRITE(*,*) "    !write!=> strss: ", reshape(virial,(/9/))

            IF (vstyle==24 .or. vstyle==25) THEN ! returns fantasy friction
                WRITE(initbuffer,'(a)') "{"
                WRITE(string, '(a,3x,f15.8,a,f15.8,a,f15.8,&
     &          3x,a)') '"dipole": [',dip(1),",",dip(2),",",dip(3),"],"
                string2 = TRIM(initbuffer) // TRIM(string)
                initbuffer = TRIM(string2)

                WRITE(string,'(a)') '"friction": ['
                string2 = TRIM(initbuffer) // TRIM(string)
                initbuffer = TRIM(string2)
                DO i=1,3*nat
                    WRITE(string,'("[ ",*(f15.8,","))') friction(i,:)
                    length = LEN_TRIM(string)
                    trimmed = TRIM(string)
                    IF(i==3*nat) THEN
                        string = TRIM(trimmed(:length-1)) // "]"
                    ELSE
                        string = TRIM(trimmed(:length-1)) // "],"
                    ENDIF
                    string2 = TRIM(initbuffer) // TRIM(string)
                    initbuffer = TRIM(string2)
                END DO
                string =  TRIM(initbuffer) // ']}'
                initbuffer = TRIM(string)
                cbuf = LEN_TRIM(initbuffer)
                CALL writebuffer(socket,cbuf) ! Writes back the fantasy friction
                IF (verbose > 1) WRITE(*,*) "!write!=> extra_length:", &
     &          cbuf
                CALL writebuffer(socket,initbuffer,cbuf)
                IF (verbose > 1) WRITE(*,*) "    !write!=> extra: ",  &
     &          initbuffer
            ELSEIF (vstyle==5 .or. vstyle==6 .or. vstyle==8) THEN ! returns the dipole through initbuffer
               WRITE(initbuffer, '(a,3x,f15.8,a,f15.8,a,f15.8, &
     &         3x,a)') '{"dipole": [',dip(1),",",dip(2),",",dip(3),"]}"
               cbuf = LEN_TRIM(initbuffer)
               CALL writebuffer(socket,cbuf) ! Writes back the molecular dipole
               IF (verbose > 1) WRITE(*,*)  &
     &         "    !write!=> extra_length: ", cbuf
               CALL writebuffer(socket,initbuffer,cbuf)
               IF (verbose > 1) WRITE(*,*) "    !write!=> extra: ", &
     &         initbuffer
            ELSEIF (vstyle==26 .OR. vstyle==3) THEN ! returns the dipole and its derivative through initbuffer
               WRITE(string, '(a,3x,f15.8,a,f15.8,a,f15.8, 3x,a)') '{"dipole": [',dip(1),",",dip(2),",",dip(3),"],"
               WRITE(string2, *) "(a,3x,", 3*nat - 1, '(f15.8, ","),f15.8,3x,a)'
               !WRITE(out_string, string2) '"dipole_derivative": [',dip_der,"]}"
               WRITE(out_string, string2) '"dipole_derivative": [',TRANSPOSE(dip_der),"],"
               out_string = TRIM(string)//TRIM(out_string)
               WRITE(string, '(a,3x, 8(f15.8, ","),f15.8,3x,a)') '"polarizability": [',polar,"]}"
               out_string = TRIM(out_string)//TRIM(string)
               cbuf = LEN_TRIM(out_string)
               CALL writebuffer(socket,cbuf) ! Writes back the molecular dipole
               CALL writebuffer(socket,TRIM(out_string),cbuf)
            ELSEIF (vstyle==1) THEN ! returns the polarizability through initbuffer
               WRITE(string, '(a,3x,f15.8,3x,a)') '{"polarizability": [',pol,"],"
               WRITE(string2, *) "(a,3x,", 3*nat - 1, '(f15.8, ","),f15.8,3x,a)'
               WRITE(out_string, string2) '"polarizability_derivative": [',pol_der,"]}"
               out_string = TRIM( TRIM(string)//TRIM(out_string) )
               cbuf = LEN_TRIM(out_string)
               CALL writebuffer(socket,cbuf) ! Writes back the gradient of total polarizability
               CALL writebuffer(socket,out_string,cbuf)
            ELSE
               cbuf = 1 ! Size of the "extras" string
               CALL writebuffer(socket,cbuf) ! This would write out the "extras" string, but in this case we only use a dummy string.
               IF (verbose > 1) WRITE(*,*)  &
     &         "    !write!=> extra_length: ", cbuf
               CALL writebuffer(socket,' ',1)
               IF (verbose > 1) WRITE(*,*)  &
     &         "    !write!=> extra: empty"
            ENDIF
            hasdata = .false.
         ELSE
            WRITE(*,*) " Unexpected header ", header
            STOP "ENDED"
         ENDIF
      ENDDO
      IF (nat > 0) DEALLOCATE(atoms, forces, msgbuffer, friction)
 
    CONTAINS
      SUBROUTINE helpmessage
         ! Help banner
         WRITE(*,*) " SYNTAX: driver.x [-u] -a address -p port -m [dummy|gas|lj|sg|harm|harm3d|morse|zundel|qtip4pf|pswater|lepsm1|lepsm2|qtip4p-efield|eckart|ch4hcbe|ljpolymer|MB|doublewell|doublewell_1D] "
         WRITE(*,*) "         -o 'comma_separated_parameters' [-v] "
         WRITE(*,*) ""
         WRITE(*,*) " For LJ potential use -o sigma,epsilon,cutoff "
         WRITE(*,*) " For SG potential use -o cutoff "
         WRITE(*,*) " For 1D/3D harmonic oscillator use -o k "
         WRITE(*,*) " For 1D morse oscillator use -o r0,D,a"
         WRITE(*,*) " For qtip4pf-efield use -o Ex,Ey,Ez with Ei in V/nm"         
         WRITE(*,*) " For ljpolymer use -o n_monomer,sigma,epsilon,cutoff "
         WRITE(*,*) " For the ideal gas, qtip4pf, zundel, ch4hcbe, nasa, doublewell or doublewell_1D no options are needed! "
       END SUBROUTINE helpmessage

   END PROGRAM
