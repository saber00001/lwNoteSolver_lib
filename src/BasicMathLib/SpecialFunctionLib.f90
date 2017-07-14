module SpecialFunctionLib
use constants
implicit none

    private
    public:: factorial
    public:: factorial_stirling
    
    
!-------------------------------------------
    interface factorial
        procedure::  factorial_1n
        procedure::  factorial_kn
    end interface factorial
    
    interface factorial_stirling
        procedure::  factorial_1n_stirling
    end interface factorial_stirling

contains

    pure function factorial_1n(n) result(factorial)
    integer(ip),intent(in)::    n
    integer(ip)::               factorial,i
        factorial = 1
        do i=1,n
            factorial = factorial*i
        enddo
    end function factorial_1n
    
    pure function factorial_kn(k,n) result(factorial)
    integer(ip),intent(in)::    k,n
    integer(ip)::               factorial,i
        factorial = 1
        do i=k,n
            factorial = factorial*i
        enddo
    end function factorial_kn
    
    !-------when n>100, fast factorial
    pure function factorial_1n_stirling(n) result(factorial)
    integer(ip),intent(in)::    n
    integer(ip)::               factorial
        factorial = nint( sqrt(2.d0*pi*n) * (n/e)**n )
    end function factorial_1n_stirling
    
end module SpecialFunctionLib