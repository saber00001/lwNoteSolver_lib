module GasChemReactionModel_
use constants
use GasPerfectThermoModel_
implicit none

    private          
    public:: GasChemReactionModel
    public:: GasSpeciesModel
    !--
    
    
    !----------------------------------------------------------------
    type GasSpeciesModel
        
        private
        real(rp),allocatable,dimension(:)::     CH
        real(rp),allocatable,dimension(:)::     CL
        real(rp)::                              mw_ ![g/mol]
            
    contains  
          
        generic::                   init    =>  init_n,     &
                                                init_coefs
        procedure::                 init_n
        procedure::                 init_coefs

        !menber function
        procedure::                 thermoFitCoefH
        procedure::                 thermoFitCoefL
        procedure::                 mw => molecularWeight
        procedure::                 MolecularWeight ![g/mol]
        procedure::                 R               ![J][K-1][g-1]
        
        procedure::                 Cp => CpFit
        procedure::                 CpFit
        
        procedure::                 H => EnthalpyFit
        procedure::                 EnthalpyFit
        
        procedure::                 S => EntropyFit
        procedure::                 EntropyFit
        
    end type GasSpeciesModel
    
    !------------------------------!
    type GasChemReactionModel
        
        private
        type(GasSpeciesModel),allocatable,dimension(:)::  species
    
    contains
    
        generic::                   init => init_sp
        procedure::                 init_sp
        
        !chemistry production of basic model, see chemKin of Ansys
        generic::                   omega => ProductionRate_Basic
        generic::                   ProductionRate => ProductionRate_Basic
        procedure,nopass::          ProductionRate_Basic
        
        !
        procedure,nopass::          ProgressRate
        
        
        !calculate rate constant => k = a T^b \exp{- \frac{Ea}{R T}}
        generic::                   k => Arrhenius
        generic::                   RateConstant => Arrhenius
        procedure,nopass::          Arrhenius
        
        !calculate reverse Rate Constant based on Equilibrium assumption
        generic::                   reverseRateConstant => reverseRateConstant_equilibrium!,
        procedure,nopass::          reverseRateConstant_equilibrium
        !procedure,nopass::          reverseRateConstant_Saha
        
        !calculate Equilibrium Constant based on experience formula
        procedure,nopass::          EquilibriumConstant
        
        procedure,nopass::          ThreeBodyConcentration
        
        procedure,nopass::          Temperature
        
    end type GasChemReactionModel
        
contains
!-------------------------------------------------------------------------
    !-------------
    pure subroutine init_n(this,ncoef,mw)
    class(GasSpeciesModel),intent(out)::    this
    integer(ip),intent(in)::                ncoef
    real(rp),intent(in)::                   mw
    
        allocate(this%CH(ncoef),this%CL(ncoef))
        this%CH = 0._rp; this%CL = 0._rp
        this%mw_ = mw
    
    end subroutine init_n
    !--
    pure subroutine init_coefs(this,Coef_H,Coef_L)
    class(GasSpeciesModel),intent(out)::    this
    real(rp),dimension(:),intent(in)::      Coef_H,Coef_L
        allocate(this%CH,source = Coef_H)
        allocate(this%CL,source = Coef_L)
    end subroutine init_coefs
    
    !--------------------!
    elemental real(rp) function thermoFitCoefH(this,i)
    class(GasSpeciesModel),intent(in)::     this
    integer(ip),intent(in)::                i
        thermoFitCoefH = this%CH(i)
    end function thermoFitCoefH
    
    !--------------------!
    elemental real(rp) function thermoFitCoefL(this,i)
    class(GasSpeciesModel),intent(in)::     this
    integer(ip),intent(in)::                i
        thermoFitCoefL = this%CL(i)
    end function thermoFitCoefL
    
    !--------------------
    elemental real(rp) function molecularWeight(this)
    class(GasSpeciesModel),intent(in)::     this
        MolecularWeight = this%mw_
    end function molecularWeight
    
    !--------------------
    elemental real(rp) function R(this)
    class(GasSpeciesModel),intent(in)::     this
        R = R_c / this%mw_
    end function R
    
!------------------------------------------------------------------------
    !specific heat capacity polynomial fitting
    elemental real(rp) function CpFit(this,T,R) result(shcp)
    class(GasSpeciesModel),intent(in):: this
    real(rp),intent(in)::               T,R
    real(rp),parameter::                temperature_boundary = 1000._rp
    
        if(t > temperature_boundary) then
            shcp = R*(this%thermoFitCoefH(1)*T**(-2) + this%thermoFitCoefH(2)*T**(-1) + this%thermoFitCoefH(3)&
                  + this%thermoFitCoefH(4)*t + this%thermoFitCoefH(5)*T**2&
                  + this%thermoFitCoefH(6)*T**3 + this%thermoFitCoefH(7)*T**4)
        else
            shcp = R*(this%thermoFitCoefL(1)*T**(-2) + this%thermoFitCoefL(2)*T**(-1)&
                  + this%thermoFitCoefL(3) + this%thermoFitCoefL(4)*t + this%thermoFitCoefL(5)*T**2&
                  + this%thermoFitCoefL(6)*T**3 + this%thermoFitCoefL(7)*T**4)
        end if
    
    end function CpFit
!------------------------------------------------------------------------
    !Specific heat enthalpy polynomial fitting
    elemental real(rp) function EnthalpyFit(this,T,R) result(enthalpy)
    class(GasSpeciesModel),intent(in):: this
    real(rp),intent(in)::               T,R
    real(rp),parameter::                temperature_boundary = 1000._rp
        
        if(t > temperature_boundary) then
            enthalpy= R*(-this%thermoFitCoefH(1)*T**(-1) + this%thermoFitCoefH(2)*log(T)&
                    + this%thermoFitCoefH(3)*t + this%thermoFitCoefH(4)*T**2/2._rp + this%thermoFitCoefH(5)*T**3/3._rp&
                    + this%thermoFitCoefH(6)*T**4/4._rp + this%thermoFitCoefH(7)*T**5/5._rp + this%thermoFitCoefH(8))
        else
            enthalpy= R*(-this%thermoFitCoefL(1)*T**(-1) + this%thermoFitCoefL(2)*log(T) + this%thermoFitCoefL(3)*t&
                    + this%thermoFitCoefL(4)*T**2/2._rp + this%thermoFitCoefL(5)*T**3/3._rp&
                    + this%thermoFitCoefL(6)*T**4/4._rp + this%thermoFitCoefL(7)*T**5/5._rp + this%thermoFitCoefL(8))
        end if
        
    end function EnthalpyFit
!------------------------------------------------------------------------
    !specific entropy polynomial fitting
    elemental real(rp) function EntropyFit(this,T,R) result(entropy)
    class(GasSpeciesModel),intent(in):: this
    real(rp),intent(in)::               T,R
    real(rp),parameter::                temperature_boundary = 1000._rp
        
        if(t > temperature_boundary) then
            entropy = R*(-this%thermoFitCoefH(1)*T**(-2)/2._rp - this%thermoFitCoefH(2)*T**(-1)&
                    + this%thermoFitCoefH(3)*log(T) + this%thermoFitCoefH(4)*t + this%thermoFitCoefH(5)*T**2/2._rp&
                    + this%thermoFitCoefH(6)*T**3/3._rp + this%thermoFitCoefH(7)*T**4/4._rp + this%thermoFitCoefH(9))
        else
            entropy = R*(-this%thermoFitCoefL(1)*T**(-2)/2._rp - this%thermoFitCoefL(2)*T**(-1)&
                    + this%thermoFitCoefL(3)*log(T) + this%thermoFitCoefL(4)*t + this%thermoFitCoefL(5)*T**2/2._rp&
                    + this%thermoFitCoefL(6)*T**3/3._rp + this%thermoFitCoefL(7)*T**4/4._rp + this%thermoFitCoefL(9))
        end if
        
    end function EntropyFit
!------------------------------------------------------------------------  
    
    
    
!-------------------------------------------------------------------------
    pure subroutine init_sp(this,sp)
    class(GasChemReactionModel),intent(out)::       this
    type(GasSpeciesModel),dimension(:),intent(in):: sp
        allocate(this%species,source=sp)
    end subroutine init_sp
    
    
    
    !this is a corrected version, details refer to wiki Arrhenius method
    elemental real(rp) function Arrhenius(a,b,Ea,R,T) result(k)
    real(rp),intent(in)::                                   a,b,Ea,R,T
        k = a * T**b * exp(-Ea/R/T)
    end function Arrhenius
    
    !production rate of concentraction[mol], Species concentraction[mol] = rho_i / MolecularWeight
    pure function ProductionRate_Basic(Species,stoichiometricF,stoichiometricB,&
                                        RateConstantF,RateConstantR,ThreeBodyConcentration) result(pr)
    real(rp),dimension(:),intent(in)::                      Species,RateConstantF,RateConstantR,ThreeBodyConcentration
    real(rp),dimension(:,:),intent(in)::                    stoichiometricF,stoichiometricB
    real(rp),dimension(size(Species))::                     pr
    integer(ip)::                                           si
        do si=1,size(pr)
            pr(si) = sum((stoichiometricB(:,si) - stoichiometricF(:,si)) * &
                            ProgressRate( Species,stoichiometricF,stoichiometricB,&
                            RateConstantF,RateConstantR &
                            ) * ThreeBodyConcentration(:) &    
                        )
        enddo
    end function ProductionRate_Basic
    
    !rate of progress variables
    pure function ProgressRate(Species,stoichiometricF,stoichiometricB,&
                                                            RateConstantF,RateConstantR) result(pr)
    real(rp),dimension(:),intent(in)::                      Species,RateConstantF,RateConstantR
    real(rp),dimension(:,:),intent(in)::                    stoichiometricF,stoichiometricB
    real(rp),dimension(size(RateConstantF))::               pr
    integer(ip)::                                           ri        
        do ri=1,size(pr)
            pr(ri) = RateConstantF(ri) * product(Species**stoichiometricF(ri,:))
            pr(ri) = pr(ri) - RateConstantR(ri) * product(Species**stoichiometricB(ri,:))
        enddo
    end function ProgressRate

    !reverseRateConstant
    pure function reverseRateConstant_equilibrium(rateConstant,EquilibriumConstant) result(rrc)
    real(rp),dimension(:),intent(in)::                      rateConstant,EquilibriumConstant
    real(rp),dimension(size(rateConstant))::                rrc
        rrc = rateConstant / EquilibriumConstant
    end function reverseRateConstant_equilibrium
    
    !stochiometric = stochiometricb - stochiometricf | (nr,ns)
    pure function EquilibriumConstant(stoichiometric,entropy,enthalpy,R,T) result(ec)
    real(rp),dimension(:),intent(in)::                      entropy,enthalpy
    real(rp),dimension(:,:),intent(in)::                    stoichiometric
    real(rp),intent(in)::                                   R,T
    real(rp),dimension(size(stoichiometric,dim=1))::        ec
    integer(ip)::                                           ri
    real(rp),parameter::                                    Patm=P_atm * 1.e-6_rp !pa -> 1 dynes/cm^2 = 10-6 bar
        do ri = 1,size(ec)
            !equilibrium constant constant
            ec(ri) = exp(sum(stoichiometric(ri,:)*entropy/R) - sum(stoichiometric(ri,:)*enthalpy/R/T))
            !equilibrium constant
            ec(ri) = ec(ri) * (Patm / R/T)**(sum(stoichiometric(ri,:)))
        enddo
    end function EquilibriumConstant
    
    !Three body concentration = sum(molality * Three body coefficient)
    pure function ThreeBodyConcentration(MC,ThreeBodyCoefficient) result(tbc)
    real(rp),dimension(:),intent(in)::                      MC
    real(rp),dimension(:,:),intent(in)::                    ThreeBodyCoefficient
    real(rp),dimension(size(ThreeBodyCoefficient,dim=1))::  tbc
    integer(ip)::                                           re,nr
        
        nr = size(ThreeBodyCoefficient,dim=1)
        do re =1,nr
            tbc(re) = sum(ThreeBodyCoefficient(re,:) * MC(:))
            !not three body reaction
            if(sum(ThreeBodyCoefficient(re,:))==0._rp) tbc(re)=1._rp
        end do
          
    end function ThreeBodyConcentration
    
    !caculate gas temperature |general Unit of Ri = Rc/mw [J][K-1][mol-1][mol][g-1] = [J][K-1][g-1]
    pure real(rp) function Temperature(T0,energy,Yi,Sp) result(T)
    real(rp),intent(in)::                                   T0,energy ! |energy = Cv*T
    real(rp),dimension(:),intent(in)::                      Yi
    type(GasSpeciesModel),dimension(:),intent(in)::         Sp
    integer(ip)::                                           counter
    real(rp)::                                              Tt,f,fp,energyScaledUI
    real(rp),dimension(size(Sp))::                          cpi,hi,Ri
    
        !scale J/kg -> J/g
        energyScaledUI = energy * 1.e-3_rp
        Ri = Sp%R()
        
        !newton iterative method| [f] is energy conservative function, and [fp] is the derived function
        counter=0; T = T0; Tt = 0._rp
        do while(counter<=15.or.abs(T-Tt) > GlobalEps * 1000._rp)
            counter = counter + 1
            Tt = T
            Cpi = Sp%Cp(Tt,Ri)
            Hi  = Sp%H(Tt,Ri)
            
            f  = sum(Yi*Hi) - sum(Yi*Ri)*Tt - energyScaledUI
            fp = sum(Yi*Cpi) - sum(Yi*Ri)
            T = Tt - f/fp
        end do
        
    end function Temperature
    
    
end module GasChemReactionModel_