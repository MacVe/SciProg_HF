module register
    implicit none
    private 
    public :: get_charge

    
    type :: atom_data
        character(2) :: name 
        integer :: charge
        integer :: s_orbital_e, p_orbital_e
    end type 

    type(atom_data), DIMENSION(10), public :: atom_matrix 
    
    contains 

    subroutine set_atom_matrix()
        integer :: i

        atom_matrix(1)%name = 'H'
        atom_matrix(2)%name = 'He'
        atom_matrix(3)%name = 'Li'
        atom_matrix(4)%name = 'Be'
        atom_matrix(5)%name = 'B'
        atom_matrix(6)%name = 'C'
        atom_matrix(7)%name = 'N'
        atom_matrix(8)%name = 'O'
        atom_matrix(9)%name = 'F'
        atom_matrix(10)%name = 'Ne'

        do i=1, size(atom_matrix)
            atom_matrix(i)%charge = i
            if (i<4) then
                atom_matrix(i)%s_orbital_e = i
                atom_matrix(i)%p_orbital_e = 0
            else if (i>4) then
                atom_matrix(i)%s_orbital_e = 4
                atom_matrix(i)%p_orbital_e = i-4
            end if
        end do     
        
        
    end subroutine

    function get_charge(atom_char) result(charge)
        CHARACTER, INTENT(IN) :: atom_char
        INTEGER :: charge, i 
        
        charge = 0

        do i=1, size(atom_matrix)
            if (atom_matrix(i)%name == atom_char) then
                atom_matrix(i)%charge = charge
            end if
        end do

        if (charge == 0) then
            print *, 'problem finding given atom', atom_char, 'in database'
        end if
    end function
    
 
end module 