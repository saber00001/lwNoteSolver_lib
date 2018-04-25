module laWrapperLib
use constants
use arrayOpsLib
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
    public:: inverseGeneralSquareMat
    public:: eigenTriDiagonal
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
    !call steqr(d,e,z,compz,info)
    ![d]    input the diagonal vector with size(n) and output the n eigenvalues in ascending order
    ![e]    input contains the off-diagonal elements of T with size(n-1)
    ![evc]  output n*n orthonomal eigenvector
    !Find all eigenvalues and eigenvectors of a tridiagonal matrix T
    pure subroutine eigenTriDiagonal(d,e,evc)
    real(rp),dimension(:),intent(inout)::                           d
    real(rp),dimension(:),intent(in)::                              e
    complex(rp),dimension(:,:),intent(out),optional::               evc
    
        if(present(evc)) then
            call steqr(d,e,evc)
        else
            call steqr(d,e)
        endif
    
    end subroutine eigenTriDiagonal
    
    
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
    
    
    
!dir$ if defined (without_imkl)
!dir$ else
    pure subroutine inverseGeneralSquareMat(a)
    real(rp),dimension(:,:),intent(inout)::     a
    integer(isp),allocatable,dimension(:)::     ipiv
        allocate(ipiv(max(1,min(size(a,dim=1),size(a,dim=2)))))
        call getrf(a,ipiv)
        call getri(a,ipiv)
    end subroutine inverseGeneralSquareMat
!dir$ end if
    
    
!dir$ if defined (without_imkl)
    !input a [(2:n,1),(1:n,2),(1:n-1,3)]
    !refer to chasing method
    !limiting: abs(a)=>(a(1,2)>a(1,3)>0),((a(i,2)>a(i,1)+a(i,3)),(a(n,2)>a(n,1))
    !a validation refer to 
    !https://wenku.baidu.com/view/a2065cb064ce0508763231126edb6f1aff0071d7.html
    !the lapackwrapper and chasing method seem no difference below...
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
        if(size(m,dim=1)/=size(m,dim=2)) call disableProgram
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
    
    
    
!-------------------------------------------------
    !gmres: suppose x is in the space x_0 + \mathcal(K)_m, then we have expression
    !x = x_0 + V_m y
    !to solve x for min{J = |b - A x|}, then transfer to
    !x_m = x0 + V_m y_m
    !y_m = argmin_y | \beta e_1 - \hat{H}_m y |_2
!-------------------------------------------------
    !solve AX=b, A is sparse matrix, full size (n*n)
    !nnz means number of non-zero
    !csr structure means compressed sparse row, consists of three parts(a(nnz),ia(n+1),ja(nnz))
    !a: the value of non-zero
    !ia: the start of a row of the sparse matrix in <nnz> dimension
    !ja: the location in this row, range (1:n)
    !maxInr range (1:n)
    pure subroutine gmres(a,ia,ja,b,maxOtr,maxInr, x)
    real(rp),dimension(:),intent(in)::      a,ia,ja,b
    integer(ip),intent(in)::                maxOtr,maxInr
    real(rp),dimension(:),intent(inout)::   x
    real(rp),parameter::                    rEps = 10._rp * GlobalEps
    real(rp),dimension(size(b))::           r
    real(rp),dimension(maxInr)::            c,s
    real(rp),dimension(maxInr+1)::          g,y
    real(rp),dimension(size(b),maxInr+1)::  v
    real(rp),dimension(maxInr+1,maxInr)::   h
    integer(ip)::                           i,j,iter,jcopy
    real(rp)::                              beta,eps,av,t,mu
    
        !outer loop
        do iter = 1,maxOtr
        
            r = b - Ax(a,ia,ja,x)
            beta = norm2(r)
            if(iter==1) eps = beta * rEps
            v(:,1) = r / beta
            
            g = 0._rp; g(1) = beta
            h = 0._rp
            
            !inner loop
            do j=1,maxInr
            
                jcopy = j
            
                v(:,j+1) = Ax(a,ia,ja,v(:,j))
                av = norm2(v(:,j+1))
                
                !--orth
                do i=1,j
                    h(i,j) = v(:,j+1) .ip. v(:,i)
                    v(:,j+1) = v(:,j+1) - h(i,j)*v(:,i)
                enddo
                h(j+1,j) = norm2(v(:,j+1))
                
                !special deal for singular sparse matrix
                if(h(j+1,j)/av < rEps) then
                    do i=1,j
                        t = v(:,j+1) .ip. v(:,i)
                        h(i,j) = h(i,j) + t
                        v(:,j+1) = v(:,j+1) - t*v(:,i)
                    enddo
                    h(j+1,j) = norm2(v(:,j+1))
                endif
                
                if(h(j+1,j)/=0._rp) v(:,j+1) = v(:,j+1) / h(j+1,j)
                
                !rot
                if(j>1) then
                    do i=1,j-1
                        h(i:i+1,j) = rot2([c(i),s(i)],h(i:i+1,j))
                    enddo
                endif
                
                mu = sqrt(sum(h(j:j+1,j)**2))
                c(j) = h(j,j)/mu
                s(j) = -h(j+1,j)/mu
                h(j,j) = c(j) * h(j,j) - s(j) * h(j+1,j)
                h(j+1,j) = 0._rp
                
                g(j:j+1) = rot2([c(j),s(j)] , g(j:j+1))
                beta = abs(g(j+1))
                
                if(beta < eps) exit
                
            enddo
            
            j = jcopy - 1
            y(j+1) = g(j+1)/h(j+1,j+1)
            
            do i=j,1,-1
                y(i) = (g(i) - (h(i,i+1:j+1) .ip. y(i+1:j+1))) / h(i,i)
            enddo
            
            do i=1,size(b)
                x(i) = x(i) +  (v(i,1:j+1) .ip. y(1:j+1))
            enddo
            
            if(beta < eps) exit
            
        enddo
        
    contains
    
        !compute Ax with csr form of A
        pure function Ax(a,ia,ja,x)
        real(rp),dimension(:),intent(in)::  a,ia,ja,x
        real(rp),dimension(size(x))::       Ax
        integer(ip)::                       i,j,k
            Ax = 0._rp
            do i=1,size(r)
                j = ia(i);
                k = ia(i+1) - 1
                Ax(i) = Ax(i) + (a(j:k) .ip. x(ja(j:k)))
            enddo
        end function Ax

        
        !--
        pure subroutine Arnoldi(v,H,ksp)
        real(rp),dimension(:),intent(in)::      v
        real(rp),dimension(:,:),intent(inout):: H,ksp
        
            
        
        end subroutine Arnoldi
        
        
        !--
        pure subroutine HouseholderArnoldi()
        !real(rp),dimension(:),intent(in)::      
        
        
        
        
        end subroutine HouseholderArnoldi
        
    end subroutine gmres

end module laWrapperLib
