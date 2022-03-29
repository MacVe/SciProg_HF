module Itterations
    implicit none

    private
    public :: give_Fmatrix

    contains

    subroutine give_Fmatrix(ao_integrals,coreH, D, n_AO, F)
        real(8), INTENT(IN)  :: ao_integrals(:,:,:,:), D(:,:), coreH(:,:)
        integer, INTENT(IN)  :: n_AO
        real(8), INTENT(OUT) :: F(:,:)
        integer :: mu, nu
        
        do nu = 1, n_AO
          do mu = 1, n_AO
             F = coreH + (((2.D0 * ao_integrals(:,:,mu,nu)) - ao_integrals(:,nu,mu,:)) * D)
          end do   
        end do
       
    end subroutine 
end 