module LAwrapperLib
use constants
!for now i use lapck by default. once i can write some simple LAsolvers, this module can be compiled conditionally
!!dir$ if defined(lapackwrapper)
use lapack95
!!dir$ end if
implicit none

    private
    public:: solvetridiagonal    
!!dir$ if defined(lapckwrapper)
    public:: solveGLES
!!dir$ end if
contains

    !getrf is short for triangular factorization of general real matrices
    !solve the general linear equations system by triangular factorization
!!dir$ if defined(lapckwrapper)
    pure subroutine solveGLES(a,b)
    real(rp),dimension(1:,1:),intent(inout)::   a
    real(rp),dimension(1:),intent(inout)::      b
    integer(isp),allocatable,dimension(:)::     ipiv
        !refer to ?getrf of lapack
        allocate(ipiv( max(1,min(size(a,dim=1),size(a,dim=2))) ))
        call getrf(a,ipiv)
        call getrs(a,ipiv,b)
    end subroutine solveGLES
!!dir$ end if
    !-----------------


    
    !input a [(2:n,1),(1:n,2),(1:n-1,3)]
    !refer to chasing method
    !limiting: abs(a)=>(a(1,2)>a(1,3)>0),((a(i,2)>a(i,1)+a(i,3)),(a(n,2)>a(n,1))
    !a validation refer to 
    !h ttps://wenku.baidu.com/view/a2065cb064ce0508763231126edb6f1aff0071d7.html
    !a wrapper of MKL or Optimization is better
    !it seems no difference below...
!dir$ if defined (lapackwrapper)
    pure subroutine solvetridiagonal(a,b)
    real(rp),dimension(1:,1:),intent(inout)::a
    real(rp),dimension(1:),intent(inout)::  b
        call dttrfb(a(2:n,1),a(1:n,2),a(1:n-1,3))
        call dttrsb(a(2:n,1),a(1:n,2),a(1:n-1,3),b)
    end subroutine solvetridiagonal
!dir$ else
    pure subroutine solvetridiagonal(a,b)
    real(rp),dimension(1:,1:),intent(in)::  a
    real(rp),dimension(1:),intent(inout)::  b
    real(rp),dimension(:),allocatable::     beta
    integer(ip)::                           i,n
        n = size(b)
        allocate(beta(n-1))
        !step1
        beta(1) = a(1,3) / a(1,2)
        do i=2,n-1
            beta(i) = a(i,3) / ( a(i,2) - a(i,1) * beta(i-1) )
        enddo
        !step2
        b(1) = b(1) / a(1,2)
        do i=2,n
            b(i) = ( b(i) - a(i,1) * b(i-1) ) / ( a(i,2) - a(i,1) * beta(i-1) )
        enddo
        !step3
        b(n) = b(n)
        do i=n-1,1,-1
            b(i) = b(i) - beta(i) * b(i+1)
        enddo
    end subroutine solvetridiagonal
!dir$ end if
    
    
end module LAwrapperLib