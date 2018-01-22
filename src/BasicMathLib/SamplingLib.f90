!This lib is used to sample from a prescribed distribution
    
module samplinglib
use constants
implicit none

    private
    
    public::    sampleUniform
    public::    sampleAcceptRejection
    public::    samplePullin                                    ! Generation of normal variates with given sample mean and variance
    public::    sampleGauss
    
    
!--------------------------------------------------
    interface sampleAcceptRejection                             ! the acceptance-rejection method
        procedure:: sampleAR_discrete
        procedure:: sampleAR_continuous
    end interface 
    
    
!--------------------------------------------------
    abstract interface
        pure real(rp) function abs_df1(i)
        import:: ip,rp
        integer(ip),intent(in)::    i
        end function abs_df1
    end interface
    
contains

!-------------------------------------------------
    pure real(rp) function sampleUniform()
!dir$ if defined(lwRandom)
    integer(ip),intent(in)::    idum
    integer(ip)::               mbig,mseed,mz,iff,ma,mj,mk,inext,inextp,i,k,ii
    real(rp)::                  fac,rf
    save ma,inext,inextp
    parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
    dimension ma(55)
    data iff/0/
    if(idum < 0.or.iff == 0) then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk<mz) mk=mk+mbig
            mj=ma(ii)
        end do
        do k=1,4
            do i=1,55
                ma(i)=ma(i)-ma(1+mod(i+30,55))
                if(ma(i)<mz) ma(i)=ma(i)+mbig
            end do
        end do
        inext=0
        inextp=31
    end if
    do while(.true.)
        inext=inext+1
        if(inext==56) inext=1
        inextp=inextp+1
        if(inextp==56) inextp=1
        mj=ma(inext)-ma(inextp)
        if(mj<mz) mj=mj+mbig
        ma(inext)=mj
        rf=mj*fac
        if(rf>1.e-8.and.rf < 0.99999999) return
    end do 
!dir$ else
    real(rp),parameter::    ranf_normalized = 2._rp**31 - 1._rp
        sampleUniform = ranf() / ranf_normalized
!dir$ end if
    end function sampleUniform
    
    !--
    pure integer(ip) function sampleAR_discrete(f1,fmax,xmin,xmax) result(x)
    procedure(abs_df1)::        f1        !the probability density function that is needed to input
    real(rp),intent(in)::       fmax
    integer(ip),intent(in)::    xmin,xmax
        do; x = xmin + int(sampleUniform()*dfloat(xmax+1-xmin),kind=ip)
            if(f1(x) / fmax > sampleUniform()) exit
        end do
    end function sampleAR_discrete
    
!-------------------------------------------------
    pure real(rp) function sampleAR_continuous(f1,fmax,xmin,xmax) result(x)
    procedure(absf1)::          f1
    real(rp),intent(in)::       fmax,xmin,xmax
        do; x = xmin + sampleUniform() * (xmax-xmin)
            if(f1(x) / fmax > sampleUniform()) exit
        end do
    end function sampleAR_continuous
    
!-------------------------------------------------
    pure function samplePullin(n,um,e) result(u)
    real(rp),intent(in)::       e,um            !e is the variance and um is the mean value of u
    integer(ip),intent(in)::    n               !the number of variates 
    real(rp),dimension(n)::     u
    real(rp)::                  mu,r,b,te,ee    
    integer(ip)::               neta,j,m,i
    real(rp),pointer,dimension(:)::    eprime,v,t
        if(n > 2) then
            if(ibits(n - 1,0,1)==0) then
                mu = 0._rp
            else if (ibits(n - 1,0,1)==1) then
                mu = 0.5_rp
            end if
            r = dfloat(n - 1) / 2._rp - mu
            neta = int(r + 2._rp * mu)  
            allocate (eprime(neta),v(n + 1),t(neta - 1)) 
            eprime = e        
            v = 0._rp
            t = 0._rp
            do j=1,neta-1
                t(j) = sampleUniform() ** (1 / (neta - j - mu))
            end do
            do j=1,neta
                do m=1,j-1
                    eprime(j) = eprime(j) * t(m)
                end do
                if(j==neta) exit
                eprime(j) = eprime(j) * (1._rp - t(j))
            end do
            do m=1,int(r)
                b = 2._rp * pi * sampleUniform()
                v(2 * m) = sqrt(2._rp * eprime(m)) * cos(b)
                v(2 * m + 1) = sqrt(2._rp * eprime(m)) * sin(b)
            end do
            if(ibits(n - 1,0,1)==1) then
                if(sampleUniform() < 0.5) then
                    v(n) = sqrt(2._rp * eprime(neta))
                else
                    v(n) = -sqrt(2._rp * eprime(neta))
                end if
            end if
            u(1) = um - sqrt(dfloat(n) - 1._rp) * v(2) / sqrt(1._rp * dfloat(n))
            do i=2,n
                u(i) = u(i - 1) + (sqrt(dfloat(n) + 2._rp - dfloat(i)) * v(i) - sqrt(dfloat(n) &
                    - 1._rp * dfloat(i)) * v(i + 1)) / sqrt(dfloat(n) + 1._rp - dfloat(i))
            end do
            deallocate (eprime,v,t)
        else
            if(sampleUniform() < 0.5) then
                ee = sqrt(e)
            else
                ee = -sqrt(e)    
            end if
            u(1) = um + ee
            u(2) = 2 * um - u(1)
        end if
    end function samplePullin
    
    !----------------------------------------
    pure real(rp) function sampleGauss()
    real(rp)::                  U1,U2
    
        U1 = sampleUniform()
        U2 = sampleUniform()
        sampleGauss = sqrt(-2._rp*log(U1))*COS(2._rp*pi*U2)
    
    end function sampleGauss
    
    
end module samplinglib
