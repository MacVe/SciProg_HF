program HartreeFock

   ! Demonstration program that can be used as a starting point
   ! Lucas Visscher, March 2022
   ! Forked by: Mac Veldhuizen, March 2022

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
     real(8)  :: E_HF, answer, prev_answer
     real(8), allocatable :: F(:,:),V(:,:),T(:,:),S(:,:), C(:,:), eps(:), D(:,:), D_old(:,:), coreH(:,:)
     integer  :: i
     real     :: converge

     ! The following large array can be eliminated when Fock matrix contruction is implemented
     real(8), allocatable :: ao_integrals (:,:,:,:)
   
     ! Definition of the molecule
     call define_molecule(molecule)

     ! Definition of the GTOs
     call define_basis(molecule, ao_basis)
     n_AO = ao_basis%nao
   
     ! Definition of the number of occupied orbitals
     n_occ = get_n_occ(molecule)

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
     
     answer = 0
     prev_answer = 0

     do i = 1, 100     
      call give_Fmatrix(ao_integrals, coreH, D, n_AO, F)
   
      call solve_genev (F,S,C,eps)    
      
      D_old = D

      do lambda = 1, n_ao
        do kappa = 1, n_ao
           D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))
        end do
      end do   

      answer = abs(sum(D - D_old))

      print *, 'Change for this convergence is: ', abs(answer - prev_answer)
      if (abs(answer - prev_answer) <= 0.0001) then !set sum to variable, if statement -> print , smaller then change 
        print *, 'System converged!'
        exit
      end if

      prev_answer = answer
     end do 

     E_HF = 2.D0 * sum(coreH*D)     
     
     do lambda = 1, n_ao
        do kappa = 1, n_ao
           E_HF = E_HF + 2.D0 *  D(kappa,lambda) * sum(D*ao_integrals(:,:,kappa,lambda))
           E_HF = E_HF - 1.D0 *  D(kappa,lambda) * sum(D*ao_integrals(:,lambda,kappa,:))
       end do
     end do
   
     print*, ''
     print*, "The Hartree-Fock energy:    ", E_HF
   end

   subroutine define_molecule(molecule)
     use molecular_structure
     use register

     type(molecular_structure_t), intent(inout) :: molecule
     real(8), ALLOCATABLE :: charge(:), coord(:,:)
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
      READ (2,*) atom_char, coord (:,i)
      print *, 'Read:', atom_char, coord (:,i)
      
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
     integer :: i

     do i=1, molecule%num_atoms
       call basis_routine_atom(molecule%charge(i), molecule%coord(:,i), ao_basis)
     end do
   end subroutine


 
