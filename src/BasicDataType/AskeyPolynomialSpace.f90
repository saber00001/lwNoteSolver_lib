!this is a space(mathematics) based on the askey polynomial basis
!a function space with different measures has different optimal basis
!it is a space due to the commutative operation (+ and *) and the measurement induced by inner product
!the operational object is the coordinates of truncated basis expressed as an array
module AskeyPolynomialSpace_
use constants
use stringOpsLib
use IntegrationLib
use laWrapperLib
use polynomial_
implicit none

    private
    public:: AskeyPolynomialSpace
    
    !--
    character(10),dimension(6),parameter:: polyType = ['legendre','hermite',    &
            'jacobi','chebyshev1','chebyshev2','laguerre']
    
!----------------------------------
    !due to the normal-orthogonal polynomial basis, assert: basis(0) = 1 => phi(0) = mean(phi)
    !details see my TexNotes/UQ/fragment/method
    type:: AskeyPolynomialSpace
    
        private
        
        integer(ip)::                               quadNp_
        integer(ip)::                               truncOd_
        
        type(polynomial),dimension(:),allocatable:: basis_
        character(10)::                             basisType_
        
        !the coordinates is calculated by the innerproduct of this space
        !the innerproduct is based on the quadrature rule
        !and they maybe differ with a coef due to the normalized measurement/probability
        real(rp)::                                  ipMeasCoef_
        
        !(0:so,1:np), use quadrature rule for inner product
        real(rp),dimension(:,:),allocatable::       ipQuadxPoly_
        real(rp),dimension(:),allocatable::         ipQuadw_
        
        !use tribasis quadrature value for multiplication and division
        !val_{ijk} = measCoef * \int \phi_i \phi_j \phi_k dp(\xi)
        real(rp),dimension(:,:,:),allocatable::     triBasisQuadVal_
        
    contains
    
        procedure::             init
        procedure::             makeOpsCoef

        !operation in space, short for use
        !--binary
        procedure::             mt
        procedure::             dv
        !--unitary
        procedure::             op
        procedure::             sqt
        procedure::             ex
        procedure::             ln
        procedure::             pw
        procedure::             dc
        
        !--derivative
        !\frac{\partial f_j}{\partial x_k} = \int \frac{\partial f}{\partial x} \frac{\partial x}{\partial x_k} 
        !\phi_j = \sum (\frac{\partial f}{\partial x})_i S_ijk
        procedure::             deri
        !diff = df = df/du du 
        procedure::             diff
        
        !member function
        procedure::             truncOd
        procedure::             quadNp
    
    end type AskeyPolynomialSpace
    
contains

    subroutine init(this,basisType,basis)
    class(AskeyPolynomialSpace),intent(out)::   this
    character(*),intent(in)::                   basisType
    type(polynomial),dimension(0:),intent(in):: basis
    integer(ip)::                               truncOrder
    character(len(basisType))::                 t
    integer(ip)::                               i
    
        !--
        truncOrder      = ubound(basis,1)
        this%truncOd_   = truncOrder
        this%quadNp_    = 4*truncOrder  !!3*so/2+1 is safe for triBasis, but not enough for other operator
        
        !--
        allocate(this%basis_(0:this%truncOd_))
        this%basis_     = basis
    
        !--
        t = basisType
        call lowerString(t)
        call check
        this%basisType_ = t
        
        select case(t)
        case(polyType(1))
            this%ipMeasCoef_ = 0.5_rp
        case(polyType(2))
            this%ipMeasCoef_ = 1._rp !to be decided
        case default
            stop 'error: AskyPolynomialSpace/init impossible error'
        end select
        
        
        !a basisType alwasy correspond to a quadrature type
        call this%makeOpsCoef(basis,basisType,this%quadNp_,this%ipMeasCoef_)
        
    contains
    
        subroutine check
            do i=1, size(polyType)
                if(t == polyType(i)) exit
            end do
            if(i==size(polyType)+1) &
            stop 'error: AskeyPolynomialSpace get an unrecognized type'
            
            if(truncOrder < 0) &
            stop 'error: AskeyPolynomialSpace get a negative order'
        end subroutine check
        
    end subroutine init
    
    !--
    subroutine makeOpsCoef(this,basis,quadType,quadNp,ipMeasCoef)
    class(AskeyPolynomialSpace),intent(inout)::     this
    type(polynomial),dimension(0:),intent(in)::     basis
    character(*),intent(in)::                       quadType
    integer(ip),intent(in)::                        quadNp
    real(rp),intent(in)::                           ipMeasCoef
    integer(ip)::                                   od,np,i,j,k
    real(rp),dimension(:),allocatable::             x
        
        if(allocated(this%ipQuadw_)) deallocate(this%ipQuadw_)
        if(allocated(this%ipQuadxPoly_)) deallocate(this%ipQuadxPoly_)
        if(allocated(this%triBasisQuadVal_)) deallocate(this%triBasisQuadVal_)
        np = quadNp
        od = ubound(basis,1)
        
        !assume the accuracy of quadrature rule is N(clenshawcurtis eg.), guarantee the tribasis computation
        if(np<3*od) stop 'error: AskyPolynomialSpace/makeOpsCoef get a bad np or order'
        
        !do not store x, but the value of polynomial in x
        allocate(x(np) , this%ipQuadw_(np))
        call quadratureRule(quadType, x, this%ipQuadw_)

        allocate(this%ipQuadxPoly_(0:od,1:np))
        do j=1,np
            do i=0,od
                this%ipQuadxPoly_(i,j) = basis(i)%funcval(x(j))
            enddo
        enddo
        
        !tribasis method is more robust for multiplication and division than quadrature method
        allocate(this%triBasisQuadVal_(0:od,0:od,0:od))
        do k=0,od
            do j=0,k
                do i=0,j
                    this%triBasisQuadVal_(i,j,k) = ipMeasCoef * &
                        sum((this%ipQuadxPoly_(i,:)*this%ipQuadxPoly_(j,:)*this%ipQuadxPoly_(k,:))*this%ipQuadw_)
                    this%triBasisQuadVal_(i,k,j) = this%triBasisQuadVal_(i,j,k)
                    this%triBasisQuadVal_(j,i,k) = this%triBasisQuadVal_(i,j,k)
                    this%triBasisQuadVal_(j,k,i) = this%triBasisQuadVal_(i,j,k)
                    this%triBasisQuadVal_(k,i,j) = this%triBasisQuadVal_(i,j,k)
                    this%triBasisQuadVal_(k,j,i) = this%triBasisQuadVal_(i,j,k)
                enddo
            enddo
        enddo
        
    end subroutine makeOpsCoef

    
    
    
!------------------------------------------
    pure function mt(this,lhs,rhs)
    class(AskeyPolynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     lhs,rhs
    real(rp),dimension(0:ubound(lhs,1))::   mt
    integer(ip)::                           i,j,k
        mt = 0._rp
        do k=0,this%truncOd_
            do j=0,this%truncOd_
                do i=0,this%truncOd_
                    mt(k) = mt(k) + lhs(i)*rhs(j)*this%triBasisQuadVal_(i,j,k)
                enddo
            end do
        enddo
    end function mt
    
    !--
    pure function dv(this,up,down)
    class(AskeyPolynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     up,down
    real(rp),dimension(0:ubound(down,1))::  dv
    real(rp),dimension(0:ubound(down,1),0:ubound(down,1)):: mat
    integer(ip)::                           i,j
        do j=0,this%truncOd_
            do i=0,this%truncOd_
                mat(i,j) = sum(down(:) * this%triBasisQuadVal_([0:this%truncOd_],j,i))
            enddo
        enddo
        dv = up; call solveGeneralLES(mat,dv)
    end function dv
    
    !--unitary operation o_k = \int func(u) \phi_k dp(\xi) = \sum func(\xi_i) \phi_k(\xi_i) w(i)
    !--if let func = df/du = df/du(u), it can be used to construct the derivative
    !--df^i/du_j = this%deri(this%op(a,df/du)) where u = \sum a_i \phi_i
    pure function op(this,a,func) result(o)
    class(AskeyPolynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     a
    procedure(absf1)::                      func
    real(rp),dimension(0:ubound(a,1))::     o
    real(rp),dimension(:),allocatable::     quadKernel
    integer(ip)::                           i,n,np
        n = this%truncOd_; np = this%quadNp_
        allocate(quadKernel(np))
        do i=1,np
            quadKernel(i) = func(sum(a(:) * this%ipQuadxPoly_(:,i))) * this%ipQuadw_(i)
        enddo
        do i=0,n
            o(i) = this%ipMeasCoef_ * sum(quadKernel * this%ipQuadxPoly_(i,:))
        enddo
    end function op
    
    !--
    pure function sqt(this,a) result(o)
    class(AskeyPolynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     a
    real(rp),dimension(0:ubound(a,1))::     o
        !solve \int \sqrt(\sum_{k=0}^{n} u_k*\psi_k) \psi_i \pi dx
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            ker = sqrt(x)
        end function ker
    end function sqt
    
    !--
    pure function pw(this,a,s) result(o)
    class(AskeyPolynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     a
    real(rp),intent(in)::                   s
    real(rp),dimension(0:ubound(a,1))::     o
        !solve \int (\sum_{k=0}^{n} u_k * \psi_k)**s \psi_i \pi dx
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            ker = x**s
        end function ker
    end function pw
    
    !--
    pure function ex(this,a) result(o)
    class(AskeyPolynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     a
    real(rp),dimension(0:ubound(a,1))::     o
        !solve \int e^(\sum_{k=0}^{n} c_k*\psi_k) \psi_i \pi dx
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            ker = exp(x)
        end function ker
    end function ex
    
    !--
    pure function ln(this,a) result(o)
    class(AskeyPolynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     a
    real(rp),dimension(0:ubound(a,1))::     o
        !solve \int log(\sum_{k=0}^{n} c_k*\psi_k) \psi_i \pi dx
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            ker = log(x)
        end function ker
    end function ln
    
    !--
    pure function dc(this,a,lowerfunc,upperfunc,dispt) result(o)
    class(AskeyPolynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     a
    procedure(absf1)::                      lowerfunc,upperfunc
    real(rp),intent(in)::                   dispt
    real(rp),dimension(0:ubound(a,1))::     o
        o = this%op(a,ker)
    contains
        elemental real(rp) function ker(x)
        real(rp),intent(in):: x
            if(x<dispt) then
                ker = lowerfunc(x)
            else
                ker = upperfunc(x)
            endif
        end function ker
    end function dc

    !<2014-Exploring emerging manycore architectures for 
    !uncertainty quantification through embedded stochastic Galerkin methods>
    !deri = df^i/du_j which is a matrix and lead to [df^i = matmul(df^i/du_j,du_j)]
    !df/du = \sum a_k \phi_k, df^i/du_j = \sum a_k <\phi_i \phi_j \phi_k>
    !for example | df^i/du_j = this%deri(this%op(a_u,df/du)) where u = \sum a_{ui} \phi_i
    pure function deri(this,a_dfdu)
    class(AskeyPolynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     a_dfdu
    real(rp),dimension(0:ubound(a_dfdu,1),0:ubound(a_dfdu,1)):: deri
    integer(ip)::                           i,j
        do j=0,this%truncOd_
            do i=0,this%truncOd_
                deri(i,j) = sum(a_dfdu*this%triBasisQuadVal_(i,j,:))
            enddo
        enddo
    end function deri
    
    !if we only concern the df^i = matmul(df^i/du_j,du_j), it is not necessary to compute df^i/du_j
    !the df/du = \sum a_i \phi_i is enough
    pure function diff(this,a_dfdu,a_du)
    class(AskeyPolynomialSpace),intent(in)::this
    real(rp),dimension(0:),intent(in)::     a_dfdu,a_du
    real(rp),dimension(0:ubound(a_dfdu,1))::diff
    integer(ip)::                           i,j
        diff = 0._rp
        do i=0,this%truncOd_
            do j=0,this%truncOd_
                diff(i) = diff(i) + sum(a_dfdu * this%triBasisQuadVal_(i,j,:)) * a_du(j)
            enddo
        enddo
    end function diff
    
!----------------------------
    pure integer(ip) function truncOd(this)
    class(AskeyPolynomialSpace),intent(in)::this
        truncOd = this%truncOd_
    end function truncOd
    
    pure integer(ip) function quadNp(this)
    class(AskeyPolynomialSpace),intent(in)::this
        quadNp = this%quadNp_
    end function quadNp
    
end module AskeyPolynomialSpace_