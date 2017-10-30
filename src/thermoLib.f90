!physical model branches into an extreme scale with various correlations.
!create an nested object to organize these models is a smarter way to refactor
!use the instatiated object always
module thermoLib
use constants
implicit none

    private
    !instantiation
    public:: GasModel
    

    !--
    type GasPhaseChemReactionModelType
    contains
        
        !chemistry production of two-body without any correlation
        generic::                   ProductionRate => ProductionRate_Basic
        procedure,nopass::          ProductionRate_Basic
        
        !
        procedure,nopass::          ProgressRate
        
        
        !calculate rate constant => k = a T^b \exp{- \frac{Ea}{R T}}
        generic::                   RateConstant => Arrhenius
        procedure,nopass::          Arrhenius
        
        !calculate reverse Rate Constant based on Equilibrium assumption
        generic::                   reverseRateConstant => reverseRateConstant_equilibrium!,
        procedure,nopass::          reverseRateConstant_equilibrium
        !procedure,nopass::          reverseRateConstant_Saha
        
        !calculate Equilibrium Constant based on experience formula
        procedure,nopass::          EquilibriumConstant
        
    end type GasPhaseChemReactionModelType
    
    
    !-----------------------------
    !assembly of all of models in gas
    !-----------------------------
    type GasModelType
        type(GasPhaseChemReactionModelType)::   ChemReactionModel
    contains
    
        !Sutherland Equation to calculate viscocity
        generic::                   Sutherland => Sutherland_custom,Sutherland_air
        procedure,nopass::          Sutherland_custom
        procedure,nopass::          Sutherland_air
        
        !-------------
        !idealGasModel
        !-------------
        !ideal Gas state Equation T(p,rho,R)
        procedure,nopass::          idealGasT

        !ideal gas isoEntropic flow
        procedure,nopass::          isoEntropicPressureRatio
        procedure,nopass::          isoEntropicTemperatureRatio

    end type GasModelType
    
    !instantiation
    type(GasModelType)::            GasModel
    
    

contains


    !--
    elemental real(rp) function Sutherland_custom(T,T0,T1,Miu0) result(miu)
    real(rp),intent(in)::   T,T0,T1,Miu0
        miu = (T0 + T1) / (T + T1)
        miu = miu * Miu0 * sqrt( (T / T0)**3 )
    end function Sutherland_custom
    
    !--
    elemental real(rp) function Sutherland_air(T) result(miu)
    real(rp),intent(in)::   T
        miu = GasModel%Sutherland(T,273.15_rp,110.4_rp,1.7161_rp)
    end function Sutherland_air
    
!-------------------------------------------------------------------------
    !this is a corrected version, details refer to wiki Arrhenius method
    elemental real(rp) function Arrhenius(a,b,Ea,R,T) result(k)
    real(rp),intent(in)::   a,b,Ea,R,T
        k = a * T**b * exp(-Ea/R/T)
    end function Arrhenius
    
    !production rate of concentraction[mol], Species concentraction[mol] = rho_i / MolecularWeight
    pure function ProductionRate_Basic(Species,stoichiometricF,stoichiometricB,&
                                            RateConstantF,RateConstantR) result(pr)
    real(rp),dimension(:),intent(in)::      Species,RateConstantF,RateConstantR
    real(rp),dimension(:,:),intent(in)::    stoichiometricF,stoichiometricB
    real(rp),dimension(size(Species))::     pr
    integer(ip)::                           si
        do si=1,size(pr)
            pr(si) = sum((stoichiometricB(:,si) - stoichiometricF(:,si)) * &
                            GasModel%ChemReactionModel%ProgressRate( &
                            Species,stoichiometricF,stoichiometricB,RateConstantF,RateConstantR &
                            ) &    
                        )
        enddo
    end function ProductionRate_Basic
    
    !rate of progress variables
    pure function ProgressRate(Species,stoichiometricF,stoichiometricB,&
                                            RateConstantF,RateConstantR) result(pr)
    real(rp),dimension(:),intent(in)::      Species,RateConstantF,RateConstantR
    real(rp),dimension(:,:),intent(in)::    stoichiometricF,stoichiometricB
    real(rp),dimension(size(RateConstantF))::pr
    integer(ip)::                           ri        
        do ri=1,size(pr)
            pr(ri) = RateConstantF(ri) * product(Species**stoichiometricF(ri,:))
            pr(ri) = pr(ri) - RateConstantR(ri) * product(Species**stoichiometricB(ri,:))
        enddo
    end function ProgressRate

    !reverseRateConstant
    pure function reverseRateConstant_equilibrium(rateConstant,EquilibriumConstant) result(rrc)
    real(rp),dimension(:),intent(in)::          rateConstant,EquilibriumConstant
    real(rp),dimension(size(rateConstant))::    rrc
        rrc = rateConstant / EquilibriumConstant
    end function reverseRateConstant_equilibrium
    
    !stochiometric = stochiometricb - stochiometricf | (nr,ns)
    pure function EquilibriumConstant(stoichiometric,entropy,enthalpy,R,T) result(ec)
    real(rp),dimension(:),intent(in)::          entropy,enthalpy
    real(rp),dimension(:,:),intent(in)::        stoichiometric
    real(rp),intent(in)::                       R,T
    real(rp),dimension(size(stoichiometric,dim=1))::  ec
    integer(ip)::                               ri
        do ri = 1,size(ec)
            !equilibrium constant constant
            ec(ri) = exp(sum(stoichiometric(ri,:)*entropy/R) - sum(stoichiometric(ri,:)*enthalpy/R/T))
            !equilibrium constant
            ec(ri) = ec(ri) * (P_atm / R/T)**(sum(stoichiometric(ri,:)))
        enddo
    end function EquilibriumConstant
    
    
    
!------------------------------------------------------------------------
    elemental real(rp) function idealGasT(rho,p,R) result(T)
    real(rp),intent(in)::   rho,p,R
        T = p / rho / R
    end function idealGasT
    
    !p0/p for ideal isoentropic flow
    elemental real(rp) function isoEntropicPressureRatio(ma,gm) result(pr)
    real(rp),intent(in)::   ma,gm
        pr = (1._rp + (gm-1._rp)/2._rp * ma**2 ) ** (gm/(gm-1._rp))
    end function isoEntropicPressureRatio
    
    !T0/T for ideal isoentropic flow
    elemental real(rp) function isoEntropicTemperatureRatio(ma,gm) result(TR)
    real(rp),intent(in)::   ma,gm
        TR = 1._rp + (gm-1._rp)/2._rp * ma**2
    end function isoEntropicTemperatureRatio    
    
end module thermoLib