module MonteCarloLib
use constants
use arrayOpsLib
use SamplingLib
use laWrapperLib
implicit none


    interface gmlmc
        procedure:: geoMultiLevelMonteCarlo
    end interface gmlmc

contains

    !alpha -> weak error ~ O(2^{-alpha*l})
    !beta -> variance ~ O(2^{-beta*l})
    !gamma -> sample cost ~ O(2^{alpha*l})
    !sumsl(1,:) = sum(Y)
    !sumsl(2,:) = sum(Y^2)
    pure subroutine geoMultiLevelMonteCarlo(fn,maxlvl,eps,nl0,alpha0,beta0,gamma,r)
    integer(ip),dimension(0:),intent(in)::  nl0 !the number of samples on lvl 0,1,2
    integer(ip),intent(in)::                maxlvl
    real(rp),intent(in)::                   eps,alpha0,beta0,gamma
    real(rp),intent(out)::                  r
    integer(ip)::                           i,l,l1,lvl,lvl1
    integer(ip),dimension(maxlvl)::         dnl,nl,nsmp  !
    real(rp),dimension(maxlvl)::            ml,vl,cl,x   !mean and variance
    real(rp)::                              alpha,beta,rem,sums(4),sumsl(4,maxlvl)
    integer(ip),dimension(:,:),allocatable::idx
    real(rp),dimension(:,:),allocatable::   A
    abstract interface
        pure function fn(lvl,N)
        import:: ip,rp
        integer(ip),intent(in)::    lvl,N
        real(rp),dimension(4)::     fn
        end function fn
    end interface
    !------------------------------------------------------
    
        alpha = max(0._rp,alpha0)
        beta = max(0._rp,beta0)
        
        lvl = 2; lvl1 = 3
        dnl = nl(:)
        nl = 0
        sumsl = 0._rp
        
        do while(sum(dnl)>0)
            
            do l = 0,lvl
                l1 = l+1
                if(dnl(l1)>0) then  !increase sampling number
                    sums = fn(l,dnl(l1))
                    nl(l1) = nl(l1) + dnl(l1)
                    sumsl(:,l1) = sumsl(:,l1) + sums(:)
                endif                
            enddo
            
            !evaluate
            ml(1:lvl1) = sumsl(1,1:lvl1) / real(nl(1:lvl1),kind=rp)
            vl(1:lvl1) = max(0._rp , sumsl(2,1:lvl1)/real(nl(1:lvl1),kind=rp) - ml(1:lvl1)**2)
            
            !fix
            do l = 3,lvl1
                ml(l) = max(ml(l) , 0.5_rp*ml(l-1)/2._rp**alpha)    !extra
                vl(l) = max(vl(l) , 0.5_rp*vl(l-1)/2._rp**beta)     !extra
            enddo
            
            !linear regressian estimate alpha and beta
            if(alpha0<=0._rp) then
                call repmat(real([1:lvl],kind=rp), 1 , 2 , A)
                call repmat(trans([1:0:-1]), lvl, 1 , idx)
                A = A**idx
                x(1:lvl) = log2(ml(2:lvl1))
                call solveGeneralLES(A,x(1:lvl))
                alpha = max(0.5_rp, -x(1))
            endif
            
            if(beta0<=0._rp) then
                call repmat(real([1:lvl],kind=rp), 1 , 2 , A)
                call repmat(trans([1:0:-1]), lvl, 1 , idx)
                A = A**idx
                x(1:lvl) = log2(vl(2:lvl1))
                call solveGeneralLES(A,x(1:lvl))
                beta = max(0.5_rp, -x(1))
            endif
            
            !--
            cl(1:lvl1) = 2._rp**(gamma*[0:lvl])
            nsmp(1:lvl1) = ceiling(2._rp*vl(1:lvl1)/cl(1:lvl1) * sum(sqrt(vl(1:lvl1)*cl(1:lvl1))) / eps**2)
            dnl(1:lvl1) = max(0 , nsmp(1:lvl1) - nl(1:lvl1))
            
            !--
            if(all(dnl < 0.01_rp*nl)==.true.) then
                rem = maxval(ml(lvl1-2:lvl1) * 2._rp**(alpha*[-2:0])) / (2._rp**alpha - 1._rp)
                if(rem>eps/sqrt(2._rp) .and. lvl1 < maxlvl) then
                    lvl = lvl1
                    lvl1 = lvl + 1
                    vl(lvl1) = vl(lvl) / 2**beta
                    nl(lvl1) = 0
                    sumsl(:,lvl1) = 0._rp
                    
                    cl(1:lvl1) = 2._rp**(gamma*[0:lvl])
                    nsmp(1:lvl1) = ceiling(2._rp*vl(1:lvl1)/cl(1:lvl1) * sum(sqrt(vl(1:lvl1)*cl(1:lvl1))) / eps**2)
                    dnl(1:lvl1) = max(0 , nsmp(1:lvl1) - nl(1:lvl1))
                endif
            endif
        enddo
        
        r = sum(sumsl(1,1:lvl1)/nl(1:lvl1))
        
    end subroutine geoMultiLevelMonteCarlo
    
    !---
    pure subroutine geoMultiLevelMonteCarloTest(fn,costfactor,N,L,N0,eps)
    real(rp),intent(in)::       costfactor,eps
    integer(ip),intent(in)::    n,l,n0
    abstract interface
        pure function fn(lvl,N)
        import:: ip,rp
        integer(ip),intent(in)::lvl,N
        real(rp),dimension(4):: fn
        end function fn
    end interface
    
        
    
    end subroutine geoMultiLevelMonteCarloTest
    
    
end module MonteCarloLib