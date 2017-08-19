module LAwrapperLib
use constants
!dir$ if .not. defined(without_imkl)
use lapack95
!dir$ end if
implicit none

    private
    
    !have emergency solver
    public:: solveTridiagonalLES
    
    !totally rely on the lapack
!dir$ if .not. defined(without_imkl)
    public:: solveGeneralLES
    public:: solveSymmetryLES
    public:: eigenSymTriDiagonal
!dir$ end if
    
    
contains

    !getrf is short for triangular factorization of general real matrices
    !solve the general linear equations system by triangular factorization
    pure subroutine solveGeneralLES(a,b)
    real(rp),dimension(1:,1:),intent(inout)::   a
    real(rp),dimension(1:),intent(inout)::      b
    integer(isp),allocatable,dimension(:)::     ipiv
        !refer to ?getrf of lapack
        allocate( ipiv( max(1,min(size(a,dim=1),size(a,dim=2))) ) )
        call getrf(a,ipiv)
        call getrs(a,ipiv,b)
    end subroutine solveGeneralLES
    
    
    !solve the symmetry linear equations system by LU factorization
    pure subroutine solveSymmetryLES(a,b)
    real(rp),dimension(1:,1:),intent(inout)::   a
    real(rp),dimension(1:),intent(inout)::      b
    integer(isp),allocatable,dimension(:)::     ipiv
        !refer to ?getrf of lapack
        allocate( ipiv( max(1,min(size(a,dim=1),size(a,dim=2))) ) )
        call sytrf(a,ipiv=ipiv)
        call sytrs(a,b,ipiv)
    end subroutine solveSymmetryLES
    
    
    ![d] is the diagonal vector with size(n)
    ![eva] input sub_diagonal in beginning n-1 position and output n eigenvalues
    ![evc] output n*n orthonomal eigenvector
    !give the eigenvalues and eigenvectors for symmetric Tridiagnoal matrix
    pure subroutine eigenSymTriDiagonal(d,eva,evc)
    real(rp),dimension(:),intent(in)::              d
    real(rp),dimension(:),intent(inout)::           eva
    real(rp),dimension(:,:),intent(out),optional::  evc
        if(present(evc)) then
            call stev(d,eva,evc)
        else
            call stev(d,eva)
        endif
    end subroutine eigenSymTriDiagonal
    
    
    !input a [(2:n,1),(1:n,2),(1:n-1,3)]
    !refer to chasing method
    !limiting: abs(a)=>(a(1,2)>a(1,3)>0),((a(i,2)>a(i,1)+a(i,3)),(a(n,2)>a(n,1))
    !a validation refer to 
    !h ttps://wenku.baidu.com/view/a2065cb064ce0508763231126edb6f1aff0071d7.html
    !the lapackwrapper and chasing method seem no difference below...
!dir$ if defined (without_imkl)
    pure subroutine solveTridiagonalLES(a,b)
    real(rp),dimension(1:,1:),intent(in)::  a
    real(rp),dimension(1:),intent(inout)::  b
    real(rp),dimension(:),allocatable::     beta
    integer(ip)::                           i,n
        n = size(a,dim=1)
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
    end subroutine solveTridiagonalLES
!dir$ else
    pure subroutine solveTridiagonalLES(a,b)
    real(rp),dimension(1:,1:),intent(inout)::a
    real(rp),dimension(1:),intent(inout)::  b
    integer(ip)::                           n
        n = size(a,dim=1)
        call dttrfb(a(2:n,1),a(1:n,2),a(1:n-1,3))
        call dttrsb(a(2:n,1),a(1:n,2),a(1:n-1,3),b)
    end subroutine solveTridiagonalLES
!dir$ end if

end module LAwrapperLib