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
        generic::                   thermoFitCoefH => thermoFitCoefH_i,thermoFitCoefH_ptr
        procedure::                 thermoFitCoefH_i
        procedure::                 thermoFitCoefH_ptr
        
        generic::                   thermoFitCoefL => thermoFitCoefL_i,thermoFitCoefL_ptr
        procedure::                 thermoFitCoefL_i
        procedure::                 thermoFitCoefL_ptr
        
        procedure::                 mw => molecularWeight
        
        procedure::                 MolecularWeight ![g/mol]
        
        procedure::                 R               ![J/(K*g)] | R_c/mw
        
        generic::                   Cp => CpFit_R,CpFit_Rc
        procedure::                 CpFit_R         !following the dimension of input R
        procedure::                 CpFit_Rc        ![J/(K*mol)], R_c as the default R
        
        !H = \int_{0}^{T} Cp d{T}
        generic::                   H => EnthalpyFit_R, EnthalpyFit_Rc
        procedure::                 EnthalpyFit_R   !following the dimension of input R, if R_species[J/(K*g)], then output [J/g]
        procedure::                 EnthalpyFit_Rc  ![J/mol], R_c as the default R
        
        !S = \int_{0}^{T} \frac{Cp d{T}}{T}
        generic::                   S => EntropyFit_R, EntropyFit_Rc
        procedure::                 EntropyFit_R    !following the dimension of input R
        procedure::                 EntropyFit_Rc   ![J/(K*mol)], R_c as the default R
        
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
    pure subroutine init_n(this,ncoef)
    class(GasSpeciesModel),intent(out)::    this
    integer(ip),intent(in)::                ncoef
        allocate(this%CH(ncoef),this%CL(ncoef))
        this%CH = 0._rp; this%CL = 0._rp
        this%mw_ = 0._rp
    end subroutine init_n
    !--
    pure subroutine init_coefs(this,Coef_H,Coef_L,mw)
    class(GasSpeciesModel),intent(out)::    this
    real(rp),dimension(:),intent(in)::      Coef_H,Coef_L
    real(rp),intent(in)::                   mw
        allocate(this%CH,source = Coef_H)
        allocate(this%CL,source = Coef_L)
        this%mw_ = mw
    end subroutine init_coefs
    
    !--------------------!
    elemental real(rp) function thermoFitCoefH_i(this,i) result(c)
    class(GasSpeciesModel),intent(in)::     this
    integer(ip),intent(in)::                i
        c = this%CH(i)
    end function thermoFitCoefH_i
    
    function thermoFitCoefH_ptr(this) result(p)
    class(GasSpeciesModel),target,intent(in)::  this
    real(rp),pointer,dimension(:)::             p
        p => this%CH
    end function thermoFitCoefH_ptr
    
    !--------------------!
    elemental real(rp) function thermoFitCoefL_i(this,i) result(c)
    class(GasSpeciesModel),intent(in)::     this
    integer(ip),intent(in)::                i
        c = this%CL(i)
    end function thermoFitCoefL_i
    
    function thermoFitCoefL_ptr(this) result(p)
    class(GasSpeciesModel),target,intent(in)::  this
    real(rp),pointer,dimension(:)::             p
        p => this%CL
    end function thermoFitCoefL_ptr
    
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
    elemental real(rp) function CpFit_R(this,T,R) result(cp)
    class(GasSpeciesModel),intent(in):: this
    real(rp),intent(in)::               T,R
    real(rp),parameter::                T_crit = 1000._rp
    real(rp),pointer,dimension(:)::     ch,cl
        
        ch => this%thermoFitCoefH()
        cl => this%thermoFitCoefL()
        
        if(t > T_crit) then
            cp = R*(ch(1)*T**(-2) + ch(2)*T**(-1) + ch(3) + ch(4)*T + ch(5)*T**2 + ch(6)*T**3 + ch(7)*T**4)
        else
            cp = R*(cl(1)*T**(-2) + cl(2)*T**(-1) + cl(3) + cl(4)*T + cl(5)*T**2 + cl(6)*T**3 + cl(7)*T**4)
        end if
    
    end function CpFit_R
    
    elemental real(rp) function CpFit_Rc(this,T) result(cp)
    class(GasSpeciesModel),intent(in):: this
    real(rp),intent(in)::               T
        cp = this%Cp(T,R_c)
    end function CpFit_Rc
    
!------------------------------------------------------------------------
    !Specific heat enthalpy polynomial fitting
    elemental real(rp) function EnthalpyFit_R(this,T,R) result(enthalpy)
    class(GasSpeciesModel),intent(in):: this
    real(rp),intent(in)::               T,R
    real(rp),parameter::                T_crit = 1000._rp
    real(rp),pointer,dimension(:)::     ch,cl
    
        ch => this%thermoFitCoefH()
        cl => this%thermofitcoefl()
        if(t > T_crit) then
            enthalpy= R*(-ch(1)*t**(-1) + ch(2)*log(t) + ch(3)*t + ch(4)*t**2/2._rp + ch(5)*t**3/3._rp &
                            + ch(6)*t**4/4._rp + ch(7)*t**5/5._rp + ch(8))
        else
            enthalpy= R*(-cl(1)*t**(-1) + cl(2)*log(t) + cl(3)*t + cl(4)*t**2/2._rp + cl(5)*t**3/3._rp &
                            + cl(6)*t**4/4._rp + cl(7)*t**5/5._rp + cl(8))
        end if
        
    end function EnthalpyFit_R
    
    !Specific heat enthalpy polynomial fitting
    elemental real(rp) function EnthalpyFit_Rc(this,T) result(enthalpy)
    class(GasSpeciesModel),intent(in):: this
    real(rp),intent(in)::               T
        enthalpy = this%H(T,R_c)
    end function EnthalpyFit_Rc
    
    
!------------------------------------------------------------------------
    !specific entropy polynomial fitting
    elemental real(rp) function EntropyFit_R(this,T,R) result(entropy)
    class(GasSpeciesModel),intent(in):: this
    real(rp),intent(in)::               T,R
    real(rp),parameter::                T_crit = 1000._rp
    real(rp),pointer,dimension(:)::     ch,cl
    
        ch => this%thermoFitCoefH()
        cl => this%thermofitcoefl()
        
        if(t > T_crit) then
            entropy = R*(-ch(1)*t**(-2)/2._rp - ch(2)*t**(-1) + ch(3)*log(t) + ch(4)*t + ch(5)*t**2/2._rp &
                            + ch(6)*t**3/3._rp + ch(7)*t**4/4._rp + ch(9))
        else
            entropy = R*(-cl(1)*t**(-2)/2._rp - cl(2)*t**(-1) + cl(3)*log(t) + cl(4)*t + cl(5)*t**2/2._rp &
                            + cl(6)*t**3/3._rp + cl(7)*t**4/4._rp + cl(9))
        end if
        
    end function EntropyFit_R
    
    !specific entropy polynomial fitting
    elemental real(rp) function EntropyFit_Rc(this,T) result(entropy)
    class(GasSpeciesModel),intent(in):: this
    real(rp),intent(in)::               T
        entropy = this%S(T,R_c)
    end function EntropyFit_Rc
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
    !the parameters at interface are all based on SI Unit
    pure function EquilibriumConstant(stoichiometric,entropy,enthalpy,T) result(ec)
    real(rp),dimension(:),intent(in)::                      entropy,enthalpy
    real(rp),dimension(:,:),intent(in)::                    stoichiometric
    real(rp),intent(in)::                                   T
    real(rp),dimension(size(stoichiometric,dim=1))::        ec
    integer(ip)::                                           ri
    real(rp),parameter::                                    Teq = P_atm / R_c * 1.e-6_rp
    !pressure should be scaled from pa -> dynes/cm2 | p_dynes = p_atm * 10
    !universe gas constant should be scaled from J/(mol K) -> ergs/(mol K) | R_ergs = R_c * 10^7
    !then Patm / R should be scalsed like p_dynes/R_ergs = P_atm / R_c * 10^(-6) | this is expression in Chemkin
        do ri = 1,size(ec)
            ec(ri) = exp(sum(stoichiometric(ri,:)*entropy/R_c) - sum(stoichiometric(ri,:)*enthalpy/R_c/T)) !Kpi
            ec(ri) = ec(ri) * (Teq/T)**(sum(stoichiometric(ri,:))) !Kci
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
        do while(counter<=15.and.abs(T-Tt) > GlobalEps * 1000._rp)
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