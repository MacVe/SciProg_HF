program HartreeFock

   ! Demonstration program that can be used as a starting point
   ! Lucas Visscher, March 2022

   use molecular_structure
   use ao_basis
   use compute_integrals
   use diagonalization
   use itterations

     implicit none

     ! Variable containing the molecular structure
     type(molecular_structure_t) :: molecule
     ! Variable containing the atomic orbital basis
     type(basis_set_info_t) :: ao_basis

     ! Variable naming as in the description of the exercise
     integer  :: n_AO, n_occ
     integer  :: kappa, lambda
     real(8)  :: E_HF
     real(8), allocatable :: F(:,:),V(:,:),T(:,:),S(:,:), C(:,:), eps(:), D(:,:), D_old(:,:), coreH(:,:)
     integer  :: i
     real     :: converge

     ! The following large array can be eliminated when Fock matrix contruction is implemented
     real(8), allocatable :: ao_integrals (:,:,:,:)
   
     ! Definition of the molecule
     call define_molecule(molecule)

     ! Definition of the GTOs
     call define_basis(ao_basis)
     n_AO = ao_basis%nao
   
     ! Definition of the number of occupied orbitals
     n_occ = 3 ! hardwired for this demonstration program, should be set via input

     ! Compute the overlap matrix
     allocate (S(n_AO,n_AO))
     call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)

     ! Compute the kinetic matrix
     allocate (T(n_AO,n_AO))
     call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)

     ! Compute the potential matrix
     allocate (V(n_AO,n_AO))
     call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)

     ! Compute the core Hamiltonian matrix (the potential is positive, we scale with -e = -1 to get to the potential energy matrix)
     allocate (F(n_AO,n_AO))
     ALLOCATE (coreH(n_AO,n_AO))
     coreH = T - V 
     F = coreH
    
     
     ! Diagonalize the Fock matrix
     allocate (C(n_AO,n_AO))
     allocate (eps(n_AO))
     call solve_genev (F,S,C,eps)
     print*, "Orbital energies for the core Hamiltonian:",eps

     ! Form the density matrix
     allocate (D(n_AO,n_AO))
     do lambda = 1, n_ao
        do kappa = 1, n_ao
           D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))
       end do
     end do

     ALLOCATE(D_old(n_AO, n_AO))
     ALLOCATE(ao_integrals(n_AO,n_AO,n_AO,n_AO))

     call generate_2int(ao_basis,ao_integrals)

     do i = 1, 100     
      call give_Fmatrix(ao_integrals, coreH, D, n_AO, F)
   
      call solve_genev (F,S,C,eps)
  
      D_old = D

      do lambda = 1, n_ao
        do kappa = 1, n_ao
           D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))
        end do
      end do
     

      if (sum((F*D_old) - (D_old*F)) /= 0) then 
        print *, 'system converged'
        exit
      end if
     end do 

     ! Cint *, Dompute the Hartree-Fock energy (this should be modified, see the notes)
     E_HF = 2.D0 * sum(coreH*D)     
     ! Compute all 2-electron integrals
     !call generate_2int (ao_basis,ao_integrals)
     do lambda = 1, n_ao
        do kappa = 1, n_ao
           E_HF = E_HF + 2.D0 *  D(kappa,lambda) * sum(D*ao_integrals(:,:,kappa,lambda))
           E_HF = E_HF - 1.D0 *  D(kappa,lambda) * sum(D*ao_integrals(:,lambda,kappa,:))
       end do
     end do
   
     print*, "The Hartree-Fock energy:    ", E_HF
   end

   subroutine define_molecule(molecule)
     ! This routine should be improved such that an arbitrary molecule can be given as input
     ! the coordinates below are for a be-he dimer oriented along the x-axis with a bond length of 2 au
     use molecular_structure
     use register

     type(molecular_structure_t), intent(inout) :: molecule
     real(8), ALLOCATABLE :: charge(:), coord(:,:)
     real :: coor_x, coor_y, coor_z
     integer :: number_of_atoms, i
     character(2) :: atom_char

     call set_atom_matrix()

     open(unit=2, file='inputMolecule.txt', action="read")
    
     READ (2,*) number_of_atoms

     ALLOCATE(charge(number_of_atoms))
     ALLOCATE(coord(3, number_of_atoms))

     coord = 0.D0

     READ(2,*) !empty lane/not used

     do i=1, number_of_atoms
      READ (2,*) atom_char, coord(1,i), coord(2,i), coord(3,i)
      print *, 'Read:', atom_char, coord(1,i), coord(2,i), coord(3,i)
      
      charge(i) = get_charge(atom_char)
      print *, 'Charge allocated: ', charge(i)
     end do
  
     call add_atoms_to_molecule(molecule, charge, coord)
   end subroutine

   subroutine define_basis(molecule, ao_basis)
     use ao_basis
     use molecular_structure
     type(basis_set_info_t), intent(inout) :: ao_basis
     type(molecular_structure_t), intent(in) :: molecule
     type(basis_func_info_t) :: gto
     real(8) :: cur_coord(3)
     integer :: i, j, number_of_orbitals, l, exp, orbital_counter

     do i=1, molecule%num_atoms
       cur_coord = molecule%coord(:,i)
       orbital_counter = molecule%charge(i)
       if (orbital_counter>4) then !has p-orbitals
         l = 1 
         j = MOD((orbital_counter - 4),2) !returns number of used/occupied p-orbitals         
         exp = 4 !WRONG?
         do while (j>0)
          call add_shell_to_basis(ao_basis,l,(cur_coord), 4.D0)
         end do 
       end if

        l = 0 !s-orbitals
        if (orbital_counter>2) then !add 2s
          exp = 4
          call add_shell_to_basis(ao_basis,l,(cur_coord), 4.D0)
        end if
       exp = 1
       call add_shell_to_basis(ao_basis,l,(cur_coord), 1.D0) !pass coordintates from molecule
     end do
   end subroutine


 
