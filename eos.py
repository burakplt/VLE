# EOS models and calculations

from math import exp, sqrt, log
from numpy import roots as np_roots
from chemsep_operation import EosInterface as dbcall

class PR76 ():
    """Peng-Robinson 1976 EOS Solver"""
    
    def fugacity_vapor(stream):
        """PENG-ROBINSON equation of state solver for vapor phase.
        :param stream: Stream object that contains chemicals.
        """
        s = stream # Set stream variable
        cs = stream.substances # Components array
        T = s.get_temperature() # get system temperature Kelvin
        P = s.get_pressure() #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = s.get_fractions()[0] #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            kappa = 0.37464+1.54226*w-0.26992*w**2 #PR kappa value
            c = 0.45724*(R**2)*(Tc**2)/Pc #PR multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #PR alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            b = (0.07780*R*Tc)/Pc 
            return b

        def calculate_kij(c1, c2, tune):
            """Calculate binary interaction parameter.
            c1, c2 is the stream components, tune: 1.2 default
            """
            Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
            Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
            k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
            return k_ij
        
        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = calculate_kij(cs[i],cs[j], 1.2) #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, B-1, A-2*B-3*B**2, B**2+2*B-A*B] # PR Z-equation
            return max(np_roots(coefficients))# Return largest root for vapor phase calculation
        
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(comp,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            amix = calculate_amix(y,T)
            bmix = calculate_bmix(y)
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            Z = calculate_Z(A,B,T)
            ak = 0 # ak sum value for inside function
            for k in range(0,len(cs)):
                ak += y[k]* (1-calculate_kij(cs[k],comp, 1.2))* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/(sqrt(8)*B)*(2*ak/amix - b/bmix)*log((Z+2.414*B)/(Z-0.414*B))
            return exp(phi)

        def h_deperture(cs):
            """Departure enthalpy with PR EOS"""
            h_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                h_dep += (-R*T**2)*(der1-der2)/0.002*y[i]
            return h_dep

        def ig_enthalpy(cs):
            enthalpy = 0
            for i in range(0,len(cs)):
                enthalpy += dbcall.ig_enthalpy(cs[i].IdealGasHeatCapacityCp, T)*y[i]
            return enthalpy/1000 #kJ/kmol
        
        def s_deperture(cs):
            """Departure entropy with PR EOS"""
            s_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                dphi = (der1-der2)/0.002
                s_dep += (-R*(T*dphi + log(calculate_phi(cs[i],T))))*y[i]
            return s_dep # J/mol.K

        def ig_entropy(cs):
            entropy = 0
            P0 = 101325 # Reference pressure in Pa
            for i in range(0,len(cs)):
                #abs_entropy = float(cs[i].AbsEntropy)
                entropy += (dbcall.ig_entropy(cs[i].IdealGasHeatCapacityCp, T) -R*1000*log(P/P0) -R*1000*log(y[i]) )*y[i]
            return entropy/1000

        def gibbs_energy():
            return (ig_enthalpy(cs)+h_deperture(cs)) - (ig_entropy(cs)+s_deperture(cs))*T 
        
        fugacity_coefficients = []
        for i in range(0,len(cs)):
            fugacity_coefficients.append( calculate_phi(cs[i],T))

        return fugacity_coefficients

    # TODO FALSE RESULT!!!!!!!!!!!! Liquid fugacity 
    def fugacity_liquid(stream):
        """PENG-ROBINSON equation of state solver for liquid phase.
        :param stream: Stream object that contains chemicals.
        """
        s = stream # Set stream variable
        cs = stream.substances # Components array
        T = s.get_temperature() # get system temperature Kelvin
        P = s.get_pressure() #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = s.get_fractions()[0] #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            kappa = 0.37464+1.54226*w-0.26992*w**2 #PR kappa value
            c = 0.45724*(R**2)*(Tc**2)/Pc #PR multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #PR alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            
            b = (0.07780*R*Tc)/Pc 
            return b

        def calculate_kij(c1, c2, tune):
            """Calculate binary interaction parameter.
            c1, c2 is the stream components, tune: 1.2 default
            """
            Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
            Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
            k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
            return k_ij

        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = calculate_kij(cs[i],cs[j], 1.2) #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, B-1, A-2*B-3*B**2, B**2+2*B-A*B] # PR Z-equation
            roots = np_roots(coefficients)
            for root in roots:
                if root > 0 and root < max(roots):
                    min_root = root         
            return min_root # Return smallest root for vapor phase calculation
        
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(comp,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            amix = calculate_amix(y,T)
            bmix = calculate_bmix(y)
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            Z = calculate_Z(A,B,T)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-calculate_kij(cs[k],comp, 1.2))* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/(sqrt(8)*B)*(2*ak/amix - b/bmix)*log((Z+2.414*B)/(Z-0.414*B))
            return exp(phi)

        def h_deperture(cs):
            """Departure enthalpy with PR EOS"""
            h_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                h_dep += (-R*T**2)*(der1-der2)/0.002*y[i]
            return h_dep

        def ig_enthalpy(cs):
            enthalpy = 0
            for i in range(0,len(cs)):
                h_vap = dbcall.dh_vaporization(cs[i].HeatOfVaporization, T, cs[i].CriticalTemperature)
                enthalpy += (dbcall.ig_enthalpy(cs[i].IdealGasHeatCapacityCp, T) - h_vap)*y[i]
            return enthalpy/1000 #kJ/kmol
        
        def s_deperture(cs):
            """Departure entropy with PR EOS"""
            s_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                dphi = (der1-der2)/0.002
                s_dep += (-R*(T*dphi + log(calculate_phi(cs[i],T))))*y[i]
            return s_dep # J/mol.K

        def ig_entropy(cs):
            entropy = 0
            P0 = 101325 # Reference pressure in Pa
            for i in range(0,len(cs)):
                #abs_entropy = float(cs[i].AbsEntropy)
                entropy += (dbcall.ig_entropy(cs[i].IdealGasHeatCapacityCp, T) -R*1000*log(P/P0) -R*1000*log(y[i]) )*y[i]
            return entropy/1000

        def gibbs_energy():
            return (ig_enthalpy(cs)+h_deperture(cs)) - (ig_entropy(cs)+s_deperture(cs))*T   

    def phi_vapor(components, temp, pressure, fractions, kij_input = None, kij_tune=None):
        """PENG-ROBINSON equation of state solver for vapor phase.
        :param components: Array that contains chemicals.
        :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
        :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
        """
        cs = components # Components array
        T = temp # get system temperature Kelvin
        P = pressure #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = fractions #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            kappa = 0.37464+1.54226*w-0.26992*w**2 #PR kappa value
            c = 0.45724*(R**2)*(Tc**2)/Pc #PR multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #PR alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            
            b = (0.07780*R*Tc)/Pc 
            return b

        kijs = {}
        
        if kij_input == None:
            def calculate_kij(c1, c2, tune):
                """Calculate binary interaction parameter.
                c1, c2 is the stream components, tune: 1.2 default
                """
                Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
                Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
                k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
                return k_ij
            
            if kij_tune != None and kij_tune <3 and kij_tune >0:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], kij_tune)
            else:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], 1.2) #Default tune 1.2
        else:
            for i in range(0,len(cs)):
                for j in range(0,len(cs)):
                    if i==j:
                        kijs[(i,j)] = 0
                    else:
                        if kij_input.get((i,j),None):
                            if abs(kij_input.get((i,j))) < 0.3:
                                kijs[(i,j)] = kij_input[(i,j)]
                            else:
                                kijs[(i,j)] = 0
                        else:
                            kijs[(i,j)] = kijs[(j,i)]
        
        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = kijs[(i,j)] #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        
        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, B-1, A-2*B-3*B**2, B**2+2*B-A*B] # PR Z-equation
            return max(np_roots(coefficients))# Return largest root for vapor phase calculation
        
        amix = calculate_amix(y,T)
        bmix = calculate_bmix(y)
        A = calculate_A(calculate_amix(y,T),T)
        B = calculate_B(calculate_bmix(y),T)
        Z = calculate_Z(A,B,T)
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(i,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            comp = cs[i]
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-kijs[(k,i)])* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/(sqrt(8)*B)*(2*ak/amix - b/bmix)*log((Z+2.414*B)/(Z-0.414*B))
            return exp(phi)
        
        fugacity_coefficients = []
        for i in range(0,len(cs)):
            fugacity_coefficients.append( calculate_phi(i,T))

        return fugacity_coefficients  

class PR78 ():
    """Peng-Robinson 1978 EOS Solver"""
    
    def fugacity_vapor(stream):
        """PENG-ROBINSON equation of state solver for vapor phase.
        :param stream: Stream object that contains chemicals.
        """
        s = stream # Set stream variable
        cs = stream.substances # Components array
        T = s.get_temperature() # get system temperature Kelvin
        P = s.get_pressure() #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = s.get_fractions()[0] #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            if w <= 491:
                kappa = 0.37464 + 1.54226*w - 0.26992*w**2 #PR kappa value
            else:
                kappa = 0.379642 + 1.48503*w - 0.164423*w**2 + 0.016666*w**3
            
            c = 0.457235*(R**2)*(Tc**2)/Pc #PR multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #PR alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            
            b = (0.077796*R*Tc)/Pc 
            return b

        def calculate_kij(c1, c2, tune):
            """Calculate binary interaction parameter.
            c1, c2 is the stream components, tune: 1.2 default
            """
            Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
            Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
            k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
            return k_ij

        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = calculate_kij(cs[i],cs[j], 1.2) #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        
        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, B-1, A-2*B-3*B**2, B**2+2*B-A*B] # PR Z-equation
            return max(np_roots(coefficients))# Return largest root for vapor phase calculation
        
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(comp,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            amix = calculate_amix(y,T)
            bmix = calculate_bmix(y)
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            Z = calculate_Z(A,B,T)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-calculate_kij(cs[k],comp, 1.2))* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/(sqrt(8)*B)*(2*ak/amix - b/bmix)*log((Z+2.414*B)/(Z-0.414*B))
            return exp(phi)

        def h_deperture(cs):
            """Departure enthalpy with PR EOS"""
            h_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                h_dep += (-R*T**2)*(der1-der2)/0.002*y[i]
            return h_dep

        def ig_enthalpy(cs):
            enthalpy = 0
            for i in range(0,len(cs)):
                enthalpy += dbcall.ig_enthalpy(cs[i].IdealGasHeatCapacityCp, T)*y[i]
            return enthalpy/1000 #kJ/kmol
        
        def s_deperture(cs):
            """Departure entropy with PR EOS"""
            s_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                dphi = (der1-der2)/0.002
                s_dep += (-R*(T*dphi + log(calculate_phi(cs[i],T))))*y[i]
            return s_dep # J/mol.K

        def ig_entropy(cs):
            entropy = 0
            P0 = 101325 # Reference pressure in Pa
            for i in range(0,len(cs)):
                #abs_entropy = float(cs[i].AbsEntropy)
                entropy += (dbcall.ig_entropy(cs[i].IdealGasHeatCapacityCp, T) -R*1000*log(P/P0) -R*1000*log(y[i]) )*y[i]
            return entropy/1000

        def gibbs_energy():
            return (ig_enthalpy(cs)+h_deperture(cs)) - (ig_entropy(cs)+s_deperture(cs))*T 
        
    # TODO FALSE RESULT!!!!!!!!!!!! Liquid fugacity 
    def fugacity_liquid(stream):
        """PENG-ROBINSON equation of state solver for liquid phase.
        :param stream: Stream object that contains chemicals.
        """
        s = stream # Set stream variable
        cs = stream.substances # Components array
        T = s.get_temperature() # get system temperature Kelvin
        P = s.get_pressure() #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = s.get_fractions()[0] #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            if w <= 491:
                kappa = 0.37464 + 1.54226*w - 0.26992*w**2 #PR kappa value
            else:
                kappa = 0.379642 + 1.48503*w - 0.164423*w**2 + 0.016666*w**3
            c = 0.45724*(R**2)*(Tc**2)/Pc #PR multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #PR alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            
            b = (0.07780*R*Tc)/Pc 
            return b

        def calculate_kij(c1, c2, tune):
            """Calculate binary interaction parameter.
            c1, c2 is the stream components, tune: 1.2 default
            """
            Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
            Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
            k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
            return k_ij

        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = calculate_kij(cs[i],cs[j], 1.2) #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, B-1, A-2*B-3*B**2, B**2+2*B-A*B] # PR Z-equation
            roots = np_roots(coefficients)
            for root in roots:
                if root > 0 and root < max(roots):
                    min_root = root         
            return min_root # Return smallest root for vapor phase calculation
        
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(comp,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            amix = calculate_amix(y,T)
            bmix = calculate_bmix(y)
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            Z = calculate_Z(A,B,T)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-calculate_kij(cs[k],comp, 1.2))* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/(sqrt(8)*B)*(2*ak/amix - b/bmix)*log((Z+2.414*B)/(Z-0.414*B))
            return exp(phi)

        def h_deperture(cs):
            """Departure enthalpy with PR EOS"""
            h_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                h_dep += (-R*T**2)*(der1-der2)/0.002*y[i]
            return h_dep

        def ig_enthalpy(cs):
            enthalpy = 0
            for i in range(0,len(cs)):
                h_vap = dbcall.dh_vaporization(cs[i].HeatOfVaporization, T, cs[i].CriticalTemperature)
                enthalpy += (dbcall.ig_enthalpy(cs[i].IdealGasHeatCapacityCp, T) - h_vap)*y[i]
            return enthalpy/1000 #kJ/kmol
        
        def s_deperture(cs):
            """Departure entropy with PR EOS"""
            s_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                dphi = (der1-der2)/0.002
                s_dep += (-R*(T*dphi + log(calculate_phi(cs[i],T))))*y[i]
            return s_dep # J/mol.K

        def ig_entropy(cs):
            entropy = 0
            P0 = 101325 # Reference pressure in Pa
            for i in range(0,len(cs)):
                #abs_entropy = float(cs[i].AbsEntropy)
                entropy += (dbcall.ig_entropy(cs[i].IdealGasHeatCapacityCp, T) -R*1000*log(P/P0) -R*1000*log(y[i]) )*y[i]
            return entropy/1000

        def gibbs_energy():
            return (ig_enthalpy(cs)+h_deperture(cs)) - (ig_entropy(cs)+s_deperture(cs))*T   

    def phi_vapor(components, temp, pressure, fractions, kij_input = None, kij_tune=None):
        """PENG-ROBINSON 78 equation of state solver for vapor phase.
        :param components: Array that contains chemicals.
        :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
        :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
        """
        cs = components # Components array
        T = temp # get system temperature Kelvin
        P = pressure #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = fractions #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            if w <= 491:
                kappa = 0.37464 + 1.54226*w - 0.26992*w**2 #PR kappa value
            else:
                kappa = 0.379642 + 1.48503*w - 0.164423*w**2 + 0.016666*w**3
            
            c = 0.457235*(R**2)*(Tc**2)/Pc #PR multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #PR alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            
            b = (0.077796*R*Tc)/Pc 
            return b

        kijs = {}
        
        if kij_input == None:
            def calculate_kij(c1, c2, tune):
                """Calculate binary interaction parameter.
                c1, c2 is the stream components, tune: 1.2 default
                """
                Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
                Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
                k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
                return k_ij
            
            if kij_tune != None and kij_tune <3 and kij_tune >0:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], kij_tune)
            else:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], 1.2) #Default tune 1.2
        else:
            for i in range(0,len(cs)):
                for j in range(0,len(cs)):
                    if kij_input.get((i,j),None):
                        if abs(kij_input.get((i,j))) < 0.3:
                            kijs[(i,j)] = kij_input[(i,j)]
                        else:
                            kijs[(i,j)] = 0
                    else:
                        kijs[(i,j)] = kijs[(j,i)]

        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = kijs[(i,j)] #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        
        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, B-1, A-2*B-3*B**2, B**2+2*B-A*B] # PR Z-equation
            return max(np_roots(coefficients))# Return largest root for vapor phase calculation
        
        amix = calculate_amix(y,T)
        bmix = calculate_bmix(y)
        A = calculate_A(calculate_amix(y,T),T)
        B = calculate_B(calculate_bmix(y),T)
        Z = calculate_Z(A,B,T)
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(i,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            comp = cs[i]
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-kijs[(k,i)])* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/(sqrt(8)*B)*(2*ak/amix - b/bmix)*log((Z+2.414*B)/(Z-0.414*B))
            return exp(phi)
            
        fugacity_coefficients = []
        for i in range(0,len(cs)):
            fugacity_coefficients.append( calculate_phi(i,T))

        return fugacity_coefficients

class RK ():
    """Redlich-Kwong EOS Solver"""
    
    def fugacity_vapor(stream):
        """RK equation of state solver for vapor phase.
        :param stream: Stream object that contains chemicals.
        """
        s = stream # Set stream variable
        cs = stream.substances # Components array
        T = s.get_temperature() # get system temperature Kelvin
        P = s.get_pressure() #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = s.get_fractions()[0] #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            a = 0.427480* (R**2) * (Tc**2.5) /Pc
            return a

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            b = (0.086640*R*Tc)/Pc 
            return b

        def calculate_kij(c1, c2, tune):
            """Calculate binary interaction parameter.
            c1, c2 is the stream components, tune: 1.2 default
            """
            Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
            Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
            k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
            return k_ij

        def calculate_amix(y):
            """a value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = 0 #calculate_kij(cs[i],cs[j], 1.2) DEFAULT kij=0
                    ai = calculate_a(cs[i]) #ai value
                    aj = calculate_a(cs[j]) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2.5) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        def calculate_Z(A,B,T):
            coefficients = [1, -1, A-B-B**2, -A*B] # PR Z-equation
            root = np_roots(coefficients)
            return max(root)# Return largest root for vapor phase calculation
        
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(comp,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            a = calculate_a(comp)
            b = calculate_b(comp)
            amix = calculate_amix(y)
            bmix = calculate_bmix(y)
            A = calculate_A(calculate_amix(y),T)
            B = calculate_B(calculate_bmix(y),T)
            Z = calculate_Z(A,B,T)
            Ai = calculate_A(a,T)
            Bi = calculate_B(b,T)
            
            phi = Bi/B*(Z-1) - log(Z-B)+ A/B*(Bi/B - 2*(Ai/A)**0.5)*log(1+B/Z)
            return exp(phi)

        def h_deperture(cs):
            """Departure enthalpy with PR EOS"""
            h_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                h_dep += (-R*T**2)*(der1-der2)/0.002*y[i]
            return h_dep

        def ig_enthalpy(cs):
            enthalpy = 0
            for i in range(0,len(cs)):
                enthalpy += dbcall.ig_enthalpy(cs[i].IdealGasHeatCapacityCp, T)*y[i]
            return enthalpy/1000 #kJ/kmol
        
        def s_deperture(cs):
            """Departure entropy with PR EOS"""
            s_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                dphi = (der1-der2)/0.002
                s_dep += (-R*(T*dphi + log(calculate_phi(cs[i],T))))*y[i]
            return s_dep # J/mol.K

        def ig_entropy(cs):
            entropy = 0
            P0 = 101325 # Reference pressure in Pa
            for i in range(0,len(cs)):
                #abs_entropy = float(cs[i].AbsEntropy)
                entropy += (dbcall.ig_entropy(cs[i].IdealGasHeatCapacityCp, T) -R*1000*log(P/P0) -R*1000*log(y[i]) )*y[i]
            return entropy/1000

        def gibbs_energy():
            return (ig_enthalpy(cs)+h_deperture(cs)) - (ig_entropy(cs)+s_deperture(cs))*T

        print(calculate_phi(cs[0],T), calculate_phi(cs[1],T))

    def phi_vapor(components, temp, pressure, fractions, kij_input = None, kij_tune=None):
        """Redlich-Kwong equation of state solver for vapor phase.
        :param components: Array that contains chemicals.
        :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
        :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
        """
        cs = components # Components array
        T = temp # get system temperature Kelvin
        P = pressure #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = fractions #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            a = 0.427480* (R**2) * (Tc**2.5) /Pc
            return a

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            b = (0.086640*R*Tc)/Pc 
            return b

        kijs = {}
    
        if kij_input == None:
            def calculate_kij(c1, c2, tune):
                """Calculate binary interaction parameter.
                c1, c2 is the stream components, tune: 1.2 default
                """
                Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
                Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
                k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
                return k_ij
            
            if kij_tune != None and kij_tune <3 and kij_tune >0:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], kij_tune)
            else:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], 1.2) #Default tune 1.2
        else:
            for i in range(0,len(cs)):
                for j in range(0,len(cs)):
                    if kij_input.get((i,j),None):
                        if abs(kij_input.get((i,j))) < 0.3:
                            kijs[(i,j)] = kij_input[(i,j)]
                        else:
                            kijs[(i,j)] = 0
                    else:
                        kijs[(i,j)] = kijs[(j,i)]

        def calculate_amix(y):
            """a value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = kijs[(i,j)]
                    ai = calculate_a(cs[i]) #ai value
                    aj = calculate_a(cs[j]) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2.5) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        def calculate_Z(A,B,T):
            coefficients = [1, -1, A-B-B**2, -A*B] # PR Z-equation
            root = np_roots(coefficients)
            return max(root)# Return largest root for vapor phase calculation
        
        amix = calculate_amix(y)
        bmix = calculate_bmix(y)
        A = calculate_A(calculate_amix(y),T)
        B = calculate_B(calculate_bmix(y),T)
        Z = calculate_Z(A,B,T)
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(i,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            comp = cs[i]
            a = calculate_a(comp)
            b = calculate_b(comp)
            Ai = calculate_A(a,T)
            Bi = calculate_B(b,T)
            
            phi = Bi/B*(Z-1) - log(Z-B)+ A/B*(Bi/B - 2*(Ai/A)**0.5)*log(1+B/Z)
            return exp(phi)

        def h_deperture(cs):
            """Departure enthalpy with PR EOS"""
            h_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                h_dep += (-R*T**2)*(der1-der2)/0.002*y[i]
            return h_dep

        def ig_enthalpy(cs):
            enthalpy = 0
            for i in range(0,len(cs)):
                enthalpy += dbcall.ig_enthalpy(cs[i].IdealGasHeatCapacityCp, T)*y[i]
            return enthalpy/1000 #kJ/kmol
        
        def s_deperture(cs):
            """Departure entropy with PR EOS"""
            s_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                dphi = (der1-der2)/0.002
                s_dep += (-R*(T*dphi + log(calculate_phi(cs[i],T))))*y[i]
            return s_dep # J/mol.K

        def ig_entropy(cs):
            entropy = 0
            P0 = 101325 # Reference pressure in Pa
            for i in range(0,len(cs)):
                #abs_entropy = float(cs[i].AbsEntropy)
                entropy += (dbcall.ig_entropy(cs[i].IdealGasHeatCapacityCp, T) -R*1000*log(P/P0) -R*1000*log(y[i]) )*y[i]
            return entropy/1000

        def gibbs_energy():
            return (ig_enthalpy(cs)+h_deperture(cs)) - (ig_entropy(cs)+s_deperture(cs))*T

        phi = []
        for i in range(len(cs)):
            phi.append(calculate_phi(i,T))
        
        return phi

class SRK ():
    """Soave-Redlich-Kwong EOS solver class"""
   
    def vapor_phase(stream):
        """SRK equation of state solver for vapor phase.
        :param stream: Stream object that contains chemicals.
        """
        s = stream # Set stream variable
        cs = stream.substances # Components array
        T = s.get_temperature() # get system temperature Kelvin
        P = s.get_pressure() #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = s.get_fractions()[0] #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            kappa = 0.48 + 1.574*w - 0.176*w**2 #SRK kappa value
            c = 0.42747*(R**2)*(Tc**2)/Pc #SRK multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #SRK alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            b = (0.08664*R*Tc)/Pc 
            return b

        def calculate_kij(c1, c2, tune):
            """Calculate binary interaction parameter.
            c1, c2 is the stream components, tune: 1.2 default
            """
            Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
            Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
            k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
            return k_ij
        
        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = calculate_kij(cs[i],cs[j], 1.2) #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        
        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, -1, A-B-B**2, -A*B] # PR Z-equation
            return max(np_roots(coefficients))# Return largest root for vapor phase calculation
        
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(comp,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            amix = calculate_amix(y,T)
            bmix = calculate_bmix(y)
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            Z = calculate_Z(A,B,T)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-calculate_kij(cs[k],comp, 1.2))* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/B*(2*ak/amix - b/bmix)*log((Z+B)/Z)
            return exp(phi)

        def h_deperture(cs):
            """Departure enthalpy with SRK EOS"""
            h_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                h_dep += (-R*T**2)*(der1-der2)/0.002*y[i]
            return h_dep

        def ig_enthalpy(cs):
            enthalpy = 0
            for i in range(0,len(cs)):
                enthalpy += dbcall.ig_enthalpy(cs[i].IdealGasHeatCapacityCp, T)*y[i]
            return enthalpy/1000 #kJ/kmol
        
        def s_deperture(cs):
            """Departure entropy with PR EOS"""
            s_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                dphi = (der1-der2)/0.002
                s_dep += (-R*(T*dphi + log(calculate_phi(cs[i],T))))*y[i]
            return s_dep # J/mol.K

        def ig_entropy(cs):
            entropy = 0
            P0 = 101325 # Reference pressure in Pa
            for i in range(0,len(cs)):
                #abs_entropy = float(cs[i].AbsEntropy)
                entropy += (dbcall.ig_entropy(cs[i].IdealGasHeatCapacityCp, T) -R*1000*log(P/P0) -R*1000*log(y[i]) )*y[i]
            return entropy/1000

        def gibbs_energy(cs):
            Pref = 101325 #Pa reference pressure
            gibbs_dep = 0
            gibbs_ig = ig_enthalpy(cs) - T*ig_entropy(cs)
            for i in range(0,len(cs)):
                gibbs_dep += R*T*log(calculate_phi(cs[i],T))*y[i] + R*T*log(P/Pref)
            return gibbs_dep + gibbs_ig
        
        def heat_capacity():
            cp = 0
            for i in range(0,len(cs)):
                cp += dbcall.general(cs[i].IdealGasHeatCapacityCp, T)*y[i]
            return cp
        
        heatCp = heat_capacity()
        enthalpy = ig_enthalpy(cs)+h_deperture(cs)
        entropy = ig_entropy(cs)+s_deperture(cs)
        gibbs  = gibbs_energy(cs)
        print("Res. Enthalpy",h_deperture(cs))
        print("IG enthapy", ig_enthalpy(cs))
        """print("Entropy", entropy)
        print("Gibbs", gibbs)
        print("Heat capacity", heatCp)"""
        #TODO bütün thermal ve fiziksel properties geri döndürülecek
       
    def liquid_phase(stream):
        """SRK equation of state solver for vapor phase.
        :param stream: Stream object that contains chemicals.
        """
        s = stream # Set stream variable
        cs = stream.substances # Components array
        T = s.get_temperature() # get system temperature Kelvin
        P = s.get_pressure() #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = s.get_fractions()[0] #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            kappa = 0.48 + 1.574*w - 0.176*w**2 #SRK kappa value
            c = 0.42747*(R**2)*(Tc**2)/Pc #SRK multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #SRK alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            b = (0.08664*R*Tc)/Pc 
            return b

        def calculate_kij(c1, c2, tune):
            """Calculate binary interaction parameter.
            c1, c2 is the stream components, tune: 1.2 default
            """
            Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
            Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
            k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
            return k_ij

        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = calculate_kij(cs[i],cs[j], 1.2) #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        def calculate_Z(T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, -1, A-B-B**2, -A*B] # PR Z-equation
            roots = np_roots(coefficients)
            for root in roots:
                if root > 0 and root < max(roots):
                    min_root = root         
            return min_root # Return smallest root for vapor phase calculation
        
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(comp,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            amix = calculate_amix(y,T)
            bmix = calculate_bmix(y)
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            Z = calculate_Z(T)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-calculate_kij(cs[k],comp, 1.2))* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/B*(2*ak/amix - b/bmix)*log((Z+B)/Z)
            return exp(phi)

        def h_deperture(cs):
            """Departure enthalpy with PR EOS"""
            h_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                h_dep += (-R*T**2)*(der1-der2)/0.002*y[i]
            return h_dep

        def ig_enthalpy(cs):
            enthalpy = 0
            for i in range(0,len(cs)):
                enthalpy += dbcall.ig_enthalpy(cs[i].IdealGasHeatCapacityCp, T)*y[i]
            return enthalpy/1000 #kJ/kmol
        
        def s_deperture(cs):
            """Departure entropy with PR EOS"""
            s_dep = 0
            for i in range(0,len(cs)):
                temp = T + 0.001
                der1 = log(calculate_phi(cs[i], temp))
                temp = T - 0.001
                der2 = log(calculate_phi(cs[i], temp))
                dphi = (der1-der2)/0.002
                s_dep += (-R*(T*dphi + log(calculate_phi(cs[i],T))))*y[i]
            return s_dep # J/mol.K

        def ig_entropy(cs):
            entropy = 0
            P0 = 101325 # Reference pressure in Pa
            for i in range(0,len(cs)):
                #abs_entropy = float(cs[i].AbsEntropy)
                entropy += (dbcall.ig_entropy(cs[i].IdealGasHeatCapacityCp, T) -R*1000*log(P/P0) -R*1000*log(y[i]) )*y[i]
            return entropy/1000

        def gibbs_energy():
            return (ig_enthalpy(cs)+h_deperture(cs)) - (ig_entropy(cs)+s_deperture(cs))*T 
        
        def heat_capacity():
            cp = 0
            for i in range(0,len(cs)):
                cp += dbcall.general(cs[i].LiquidHeatCapacityCp, T)*y[i]
            return cp
        
        heatCp = heat_capacity()
        enthalpy = ig_enthalpy(cs) + h_deperture(cs)
        entropy = ig_entropy(cs) + s_deperture(cs)
        print("Enthalpy", enthalpy)
        #print("Entropy", entropy)
        #print("Heat Capacity", heatCp)
    
    def phi_vapor(components, temp, pressure, fractions, kij_input = None, kij_tune=None):
        """Soave-Redlich-Kwong equation of state solver for vapor phase.
        :param components: Array that contains chemicals.
        :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
        :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
        """
        cs = components # Components array
        T = temp # get system temperature Kelvin
        P = pressure #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = fractions #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            kappa = 0.48 + 1.574*w - 0.176*w**2 #SRK kappa value
            c = 0.42747*(R**2)*(Tc**2)/Pc #SRK multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #SRK alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            b = (0.08664*R*Tc)/Pc 
            return b

        kijs = {}
        
        if kij_input == None:
            def calculate_kij(c1, c2, tune):
                """Calculate binary interaction parameter.
                c1, c2 is the stream components, tune: 1.2 default
                """
                Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
                Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
                k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
                return k_ij
            
            if kij_tune != None and kij_tune <3 and kij_tune >0:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], kij_tune)
            else:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], 1.2) #Default tune 1.2
        else:
            for i in range(0,len(cs)):
                for j in range(0,len(cs)):
                    if kij_input.get((i,j),None):
                        if abs(kij_input.get((i,j))) < 0.3:
                            kijs[(i,j)] = kij_input[(i,j)]
                        else:
                            kijs[(i,j)] = 0
                    else:
                        kijs[(i,j)] = kijs[(j,i)]
        
        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = kijs[(i,j)] #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        
        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, -1, A-B-B**2, -A*B] # PR Z-equation
            return max(np_roots(coefficients))# Return largest root for vapor phase calculation
        
        amix = calculate_amix(y,T)
        bmix = calculate_bmix(y)
        A = calculate_A(calculate_amix(y,T),T)
        B = calculate_B(calculate_bmix(y),T)
        Z = calculate_Z(A,B,T)
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(i,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            comp = cs[i]
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-kijs[(k,i)])* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
            
            phi = b*(Z-1)/bmix - log(Z-B) - A/B*(2*ak/amix - b/bmix)*log((Z+B)/Z)
            return exp(phi)
        fug_phi = []
        for i in range(0,len(cs)):
            fug_phi.append( calculate_phi(i,T) )
        return fug_phi

    def phi_liquid(components, temp, pressure, fractions, kij_input = None, kij_tune=None):
        """Soave-Redlich-Kwong equation of state solver for vapor phase.
        :param components: Array that contains chemicals.
        :param kij_input: Dict object {(i,j):kij, (i,k):kik....}
        :param kij_tune: Tuning parameter for kij equation. Leave as None if kij_input given.
        """
        cs = components # Components array
        T = temp # get system temperature Kelvin
        P = pressure #get system pressure Pascal
        R = 8.314462 #Universal gas constant J/mol.K
        y = fractions #Molar fractions array
        
        #Calculate a(T) and b for each pure substance
        def calculate_a(component,T):
            """Input a substance i.e cs[i]
            Returns a value a = Pa.m^6/mol^2 """
            w = float(component.AcentricityFactor) #acentric factor
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            Tr = T/Tc #Reduced Temperature T is the global Temp value
            kappa = 0.48 + 1.574*w - 0.176*w**2 #SRK kappa value
            c = 0.42747*(R**2)*(Tc**2)/Pc #SRK multiply factor
            alfaT = (1 + kappa*(1-Tr**0.5))**2 #SRK alfa(T) function
            aT = c*alfaT # a(T) Equation
            return aT

        def calculate_b(component):
            """Input a substance cs[i]
            Returns b value b = m^3/mol """
            Tc = float(component.CriticalTemperature)
            Pc = float(component.CriticalPressure)
            b = (0.08664*R*Tc)/Pc 
            return b

        kijs = {}
        
        if kij_input == None:
            def calculate_kij(c1, c2, tune):
                """Calculate binary interaction parameter.
                c1, c2 is the stream components, tune: 1.2 default
                """
                Vc1 = float(c1.CriticalVolume) #Critical volume for substance 1
                Vc2 = float(c2.CriticalVolume) #Critical volume for substance 2
                k_ij = 1 - ( 2*sqrt( (Vc1**0.333)*(Vc2**0.333) )/(Vc1**0.333 + Vc2**0.333))**tune
                return k_ij
            
            if kij_tune != None and kij_tune <3 and kij_tune >0:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], kij_tune)
            else:
                for i in range(0,len(cs)):
                    for j in range(0,len(cs)):
                        kijs[(i,j)] = calculate_kij(cs[i],cs[j], 1.2) #Default tune 1.2
        else:
            for i in range(0,len(cs)):
                for j in range(0,len(cs)):
                    if kij_input.get((i,j),None):
                        if abs(kij_input.get((i,j))) < 0.3:
                            kijs[(i,j)] = kij_input[(i,j)]
                        else:
                            kijs[(i,j)] = 0
                    else:
                        kijs[(i,j)] = kijs[(j,i)]
        
        def calculate_amix(y,T):
            """a(T) value for mixture"""
            amix = 0 #Placeholder for a_mixture values
            
            for i in range(0,len(cs)) :
                for j in range(0,len(cs)):
                    kij = kijs[(i,j)] #kij value calculation
                    ai = calculate_a(cs[i],T) #ai value
                    aj = calculate_a(cs[j],T) #aj value
                    amix += y[i]*y[j]*sqrt(ai * aj)*(1-kij) #Update a_mix
            return amix
        
        def calculate_bmix(y):
            """ b value for the mixture"""
            bmix = 0
            for i in range(0, len(cs)):
                bmix += y[i]*calculate_b(cs[i])
            return bmix
        
        #amix = calculate_amix(y) # amix calculated value
        #bmix = calculate_bmix(y) #bmix calculated value

        def calculate_A(a,T):
            """Calculates A value for component or mixture. a or amix"""
            A = a * P/(R**2)/(T**2) # A factor
            return A
        
        def calculate_B(b,T):
            """Calculates B value for a component or mixture."""
            B = b * P/(R*T) # B factor
            return B

        def calculate_Z(A,B,T):
            A = calculate_A(calculate_amix(y,T),T)
            B = calculate_B(calculate_bmix(y),T)
            coefficients = [1, -1, A-B-B**2, -A*B] # PR Z-equation
            roots = np_roots(coefficients)
            for root in roots:
                if root > 0 and root < max(roots):
                    min_root = root         
            return min_root # Return smallest root for vapor phase calculation
        
        amix = calculate_amix(y,T)
        bmix = calculate_bmix(y)
        A = calculate_A(calculate_amix(y,T),T)
        B = calculate_B(calculate_bmix(y),T)
        Z = calculate_Z(A,B,T)
        # CALCULATE FUGACITY COEFFICIENT
        #Z = calculate_Z(A,B)
        def calculate_phi(i,T):
            """Vapor phase fugacity coefficient phi for a component.
            :param comp: Input the substance/chemical"""
            comp = cs[i]
            a = calculate_a(comp,T)
            b = calculate_b(comp)
            ak = 0 # ak sum value for inside function
            
            for k in range(0,len(cs)):
                ak += y[k]* (1-kijs[(k,i)])* sqrt(calculate_a(cs[k],T)*calculate_a(comp,T))
                
            phi = b*(Z-1)/bmix - log(Z-B) - A/B*(2*ak/amix - b/bmix)*log((Z+B)/Z)
            return exp(phi)
        fug_phi = []
        for i in range(0,len(cs)):
            fug_phi.append( calculate_phi(i,T) )
        return fug_phi

class Ideal():
    """Ideal property method"""
    def phi_vapor(components, temp, pressure, fractions, kij_input = None, kij_tune=None):

        phi = []
        for i in range(0, len(components)):
            phi.append(1)
        return phi

class LeeKesler():
    """Lee-Kesler methods"""
    
    def heat_properties(component, T, P):
        """T: temperature in K, P pressure in Pa
        Returns enthalpy, entropy, heat capacity"""
        #Parameters s: simple fluid, r: reference fluid
        b1_s = 0.1181193
        b1_r = 0.2026579
        b2_s = 0.265728
        b2_r = 0.331511
        b3_s = 0.154790
        b3_r = 0.027655
        b4_s = 0.030323
        b4_r = 0.203488
        c1_s = 0.0236744
        c1_r = 0.0313385
        c2_s = 0.0186984
        c2_r = 0.0503618
        c3_s = 0
        c3_r = 0.016901
        c4_s = 0.042724
        c4_r = 0.041577
        d1_s = 0.155488e-4
        d1_r = 0.48736e-4
        d2_s = 0.623689e-4
        d2_r = 0.0740336e-4
        beta_s = 0.65392
        beta_r = 1.226
        gama_s = 0.060167
        gama_r = 0.03754
        w_ref = 0.3978
        R = 8.3145

        Tc = float(component.CriticalTemperature)
        Pc = float(component.CriticalPressure)
        w = float(component.AcentricityFactor)
        Tr = T/Tc
        Pr = P/Pc

        def Vr_simple():
            B = b1_s - b2_s/Tr - b3_s/Tr**2 - b4_s/Tr**3
            C = c1_s - c2_s/Tr + c3_s/Tr**3
            D = d1_s + d2_s/Tr
            Vr = 0.5*Tr/Pr # Initial guess for Vr value. Assume Z = 0.5
            eps = 1e-4 # Tolerance
            counter = 0

            while True:
                def f(Vr):
                    func = 1 + B/Vr + C/Vr**2 + D/Vr**5 + c4_s/(Tr**3*(Vr**2))*(beta_s + gama_s/Vr**2)*exp(-gama_s/Vr**2) - Pr*Vr/Tr
                    return func
                
                def df(Vr):
                    dfunc = (f(Vr+0.001) - f(Vr-0.001))/0.002
                    return dfunc
                
                Vr_new = Vr - f(Vr)/df(Vr)
                if abs((Vr_new-Vr)/Vr) < eps:
                    break
                counter += 1
                Vr = Vr_new
                if counter > 50:
                    raise Exception("Vr simple couldn't converged in 50 loops!")
                    break
            return Vr_new

        def Vr_ref():
            B = b1_r - b2_r/Tr - b3_r/Tr**2 - b4_r/Tr**3
            C = c1_r - c2_r/Tr + c3_r/Tr**3
            D = d1_r + d2_r/Tr
            Vr = 0.5*Tr/Pr # Initial guess for Vr value. Assume Z = 0.5
            eps = 1e-4 # Tolerance
            counter = 0

            while True:
                def f(Vr):
                    func = 1 + B/Vr + C/Vr**2 + D/Vr**5 + c4_r/(Tr**3*(Vr**2))*(beta_r + gama_r/Vr**2)*exp(-gama_r/Vr**2) - Pr*Vr/Tr
                    return func
                
                def df(Vr):
                    dfunc = (f(Vr+0.001) - f(Vr-0.001))/0.002
                    return dfunc
                
                Vr_new = Vr - f(Vr)/df(Vr)
                if abs((Vr_new-Vr)/Vr) < eps:
                    break
                counter += 1
                Vr = Vr_new
                if counter > 100:
                    raise Exception("Vr reference couldn't converged in 100 loops!")
                    break
            return Vr_new    
        Vr_s = Vr_simple()
        Vr_r = Vr_ref()

        Z_simple = Pr*Vr_s/Tr
        Z_reference = Pr*Vr_r/Tr
        E_simple = c4_s/(2*Tr**3*gama_s)*( beta_s +1 -exp(-gama_s/Vr_s**2)*(beta_s+1+gama_s/Vr_s**2) )
        E_ref = c4_r/(2*Tr**3*gama_r)*( beta_r +1 -exp(-gama_r/Vr_r**2)*(beta_r+1+gama_r/Vr_r**2) )

        H_simple = Tr*( Z_simple-1 -(b2_s+2*b3_s/Tr+3*b4_s/Tr**2)/Tr/Vr_s- (c2_s-3*c3_s/Tr**2)/2/Tr/Vr_s**2+ d2_s/5/Tr/Vr_s**2 +3*E_simple )
        H_ref = Tr*( Z_reference-1 -(b2_r+2*b3_r/Tr+3*b4_r/Tr**2)/Tr/Vr_r- (c2_r-3*c3_r/Tr**2)/2/Tr/Vr_r**2+ d2_r/5/Tr/Vr_r**2 +3*E_ref )
        H_departure = H_simple + w/w_ref *(H_ref - H_simple)
        
        S_simple = log(Z_simple)- (b2_s+b3_s/Tr**2+2*b4_s/Tr**3)/Vr_s -(c1_s-2*b4_s/Tr**3)/2/Vr_s**2+ d1_s/5/Vr_s**5 +2*E_simple -log(P/101325)
        S_ref = log(Z_reference)- (b2_r+b3_r/Tr**2+2*b4_r/Tr**3)/Vr_r -(c1_r-2*b4_r/Tr**3)/2/Vr_r**2+ d1_r/5/Vr_r**5 +2*E_ref -log(P/101325)
        S_departure = S_simple + w/w_ref *(S_ref - S_simple)
        
        Cv_simple = 2*(b3_s+3*b4_s/Tr)/Vr_s*Tr**2- 3*c3_s/Tr**3/Vr_s**2 -6*E_simple
        Cv_ref = 2*(b3_r+3*b4_r/Tr)/Vr_r*Tr**2- 3*c3_r/Tr**3/Vr_r**2 -6*E_ref
        Cv_departure = Cv_simple + w/w_ref*(Cv_ref - Cv_simple)
        Cp_departure = Cv_departure -1

        return H_departure*R*Tc, S_departure*R, Cv_departure*R, Cp_departure*R
