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
    public:: fft,ifft
    
    
!------------------------------
    !the input should be equally-spaced samples of a function
    !compute X_k = sum_{n=1}^{N} x_n \cdot e^{- 2 \pi (k-1) (n-1) / N}, \quad k \in [1,N]
    interface fft
        procedure:: fft_1d_cdp
        procedure:: fft_1d_rdp
    end interface
    
    !compute x_k = \frac{1}{N} sum_{n=1}^{N} X_n \cdot e^{2 \pi (k-1) (n-1) / N}, \quad k \in [1,N]
    interface ifft
        procedure:: ifft_1d_rdp
        procedure:: ifft_1d_cdp
    end interface ifft
    
contains


    subroutine fft_1d_cdp(x)
    complex(rdp),dimension(:),intent(inout)::   x
    type(dfti_descriptor),pointer::             hand
    integer::                                   status
        status = dftiCreateDescriptor(hand, dfti_double, dfti_complex, 1, size(x))
        status = dftiCommitDescriptor(hand)
        status = dftiComputeForward(hand,x)
        status = dftiFreeDescriptor(hand)
    end subroutine fft_1d_cdp
    !--
    !subroutine fft_2d_cdp(x)
    !complex(rdp),dimension(:,:),intent(inout):: x
    !complex(rdp),dimension(size(x,dim=1)*size(x,dim=2))::   x1
    !type(dfti_descriptor),pointer::             hand
    !integer::                                   status,l(2)
    !equivalence(x,x1)
    !    l(1) = size(x,dim=1)
    !    l(2) = size(x,dim=2)
    !    status = dftiCreateDescriptor(hand, dfti_double, dfti_complex, 2, l)
    !    status = dftiCommitDescriptor(hand)
    !    status = dftiComputeForward(hand,x1)
    !    status = dftiFreeDescriptor(hand)
    !end subroutine fft_2d_cdp

    
    !-----------
    subroutine fft_1d_rdp(x)
    real(rdp),dimension(:),intent(inout)::      x
    complex(rdp),dimension(size(x))::           cx
    type(dfti_descriptor),pointer::             hand
    integer::                                   status
        cx = cmplx(x)
        call fft(cx)
        x = real(cx)
    end subroutine fft_1d_rdp
    
    
!-------------------------------------------
    subroutine ifft_1d_cdp(x)
    complex(rdp),dimension(:),intent(inout)::   x
    type(dfti_descriptor),pointer::             hand
    integer::                                   status,n
        n = size(x)
        status = dftiCreateDescriptor(hand, dfti_double, dfti_complex, 1, n)
        status = dftiCommitDescriptor(hand)
        status = dftiComputeBackward(hand,x)
        status = dftiFreeDescriptor(hand)
        x = x / dfloat(n)
    end subroutine ifft_1d_cdp
    

    subroutine ifft_1d_rdp(x)
    real(rdp),dimension(:),intent(inout)::      x
    complex(rdp),dimension(size(x))::           cx
    type(dfti_descriptor),pointer::             hand
    integer::                                   status,n
        cx = cmplx(x)
        call ifft(cx)
        x = real(cx) / dfloat(size(x))
    end subroutine ifft_1d_rdp
    

end module fftWrapperLib