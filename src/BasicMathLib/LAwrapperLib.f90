module LAwrapperLib
use constants
!dir$ if .not. defined(without_imkl)
use lapack95
!dir$ end if
implicit none

    private
    
    !have emergency solver
    public:: solveGeneralLES
    public:: solveTridiagonalLES
    public:: triFactorSquareMat
    
    !totally rely on the lapack
!dir$ if .not. defined(without_imkl)
    public:: solveSymmetryLES
    public:: eigenSymTriDiagonal
!dir$ end if
    
    
contains

    
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
    
    
!dir$ if defined (without_imkl)
    !because triFacotrSquareMat may lead U be singluar, this procedure may fails
    pure subroutine solveGeneralLES(a,b)
    real(rp),dimension(:,:),intent(inout)::     a
    real(rp),dimension(:),intent(inout)::       b
    real(rp),dimension(size(b))::               y
    integer(ip)::                               n,i
        call triFactorSquareMat(a)
        n = size(b)
        !solve LY=b
        y(1) = b(1)
        do i=2,n
            y(i) = b(i) - sum( a(i,1:i-1) * y(1:i-1) )
        enddo
        !solve UX=Y where b as X
        b(n) = y(n) / a(n,n)
        do i=n-1,1,-1
            b(i) = ( y(i) - sum( a(i,i+1:n) * b(i+1:n) ) ) / a(i,i)
        enddo
    end subroutine solveGeneralLES
!dir$ else
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
!dir$ end if
    
    
    !input a [(2:n,1),(1:n,2),(1:n-1,3)]
    !refer to chasing method
    !limiting: abs(a)=>(a(1,2)>a(1,3)>0),((a(i,2)>a(i,1)+a(i,3)),(a(n,2)>a(n,1))
    !a validation refer to 
    !https://wenku.baidu.com/view/a2065cb064ce0508763231126edb6f1aff0071d7.html
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
    
    
    
    
!--------------------------------------------------------------
!emergency only and very limited
    
    !m[in] is matrix, and m[out] is [L-I+U] witout the unit dialog of L
    !this is only for square matrix, and without ipiv unlike getrf in lapack
    !--
    !this simple procedure may fail, details can refer to <LU decompostion> in wiki
    !considering a11=0, then u11=0, U is singular. but A may be not singular, so a logic error 
    !a Partial Pivoting is needed for a proper implementation, and waiting for improving
    pure subroutine triFactorSquareMat(m)
    real(rp),dimension(:,:),intent(inout):: m
    integer(ip)::                           i,j,n
!dir$ if defined (lwcheck)
        if(size(m,dim=1)/=size(m,dim=2)) call disableprogram
!dir$ end if

        n = size(m,dim=1)
        m(2:n,1) = m(2:n,1) / m(1,1)
        do j=2,n-1
            !upper part | if(i==1) m(i,j) = m(i,j)
            do i=2,j
                m(i,j) = m(i,j) - sum(m(i,1:i-1)*m(1:i-1,j))
            enddo
            !lower part
            do i=j+1,n
                m(i,j) = ( m(i,j) - sum(m(i,1:j-1)*m(1:j-1,j)) ) / m(j,j)
            enddo
        enddo
        do i=2,n
            m(i,n) = m(i,n) - sum(m(i,1:i-1)*m(1:i-1,n))
        enddo
    
    end subroutine triFactorSquareMat

end module LAwrapperLib