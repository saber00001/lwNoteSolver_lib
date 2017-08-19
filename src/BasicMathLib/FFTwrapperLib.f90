!dir$ if .not. defined(without_imkl)
include 'mkl_dfti.f90'
!dir$ end if

!-------------------------------------------------------------
module fftWrapperLib
use constants
!dir$ if .not. defined(without_imkl)
use mkl_dfti
!dir$ end if
implicit none


    private
    !
    public:: ifft
    
    
!------------------------------
    interface ifft
        procedure:: ifft_1d_rdp
    end interface ifft
    
contains

    pure subroutine ifft_1d_rdp(x)
    real(rdp),dimension(:),intent(inout):: x
    
    end subroutine ifft_1d_rdp

end module fftWrapperLib