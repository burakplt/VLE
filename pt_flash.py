import eos, activity
from chemsep_operation import get_chemical
from chemsep_operation import EosInterface as T_prop
from math import exp, log


def bubble_T(components, z, P, model=None, kij_input = None, kij_tune=None):
    """ model = SRK, Ideal, NRTL, UNIQUAC, UNIFAC, MODFAC """
    c = components
    count = 0
    temp = 320 # Initial guess for temperature
    liquid_method = {"SRK":eos.SRK.phi_liquid,"Ideal":activity.Ideal.gamma, "NRTL":activity.NRTL.gamma,"Uniquac":activity.Uniquac.gamma,"Unifac":activity.Unifac.gamma,"Dortmund":activity.Dortmund.gamma}
    
    if model == None or model == "Ideal":
        while True:
            count+=1
            def f(temp):
                P_sum = 0
                for i in range(0,len(c)):
                    P_sum += T_prop.general( c[i].VaporPressure, temp )*z[i]
                return P - P_sum
            
            def df():
                return (f(temp + 0.001) - f(temp -0.001))/0.002
            
            f_temp = f(temp)
            df_temp = df()
            
            if abs(f_temp) < 0.00001:
                break
            else:
                temp = temp - f_temp/df_temp
            if count > 100:
                raise Warning("Bubble point calculation could not converged in 100 iterations!")
        return temp

    elif model == "SRK":
        while True:
            count+=1
            def f(temp):
                P_sum = 0
                phi = liquid_method[model](c,temp,P,z, kij_input, kij_tune)
                for i in range(0,len(c)):
                    P_sum += T_prop.general( c[i].VaporPressure, temp )*z[i]*phi[i]
                return P - P_sum
            
            def df():
                return (f(temp + 0.001) - f(temp -0.001))/0.002
            
            f_temp = f(temp)
            df_temp = df()
            
            if abs(f_temp) < 0.00001:
                break
            else:
                temp = temp - f_temp/df_temp
            if count > 100:
                raise Warning("Bubble point calculation not converged in 100 iterations!")
        return temp

    else:
        while True:
            count+=1
            def f(temp):
                P_sum = 0
                gamma = liquid_method[model](c,temp,z)
                for i in range(0,len(c)):
                    P_sum += T_prop.general( c[i].VaporPressure, temp )*z[i]*gamma[i]
                return P - P_sum
            
            def df():
                return (f(temp + 0.001) - f(temp -0.001))/0.002
            
            f_temp = f(temp)
            df_temp = df()
            
            if abs(f_temp) < 0.00001:
                break
            else:
                temp = temp - f_temp/df_temp
            if count > 100:
                raise Warning("Bubble point calculation not converged in 100 iterations!")
        return temp

def dew_T(components, z, P, model=None, kij_input = None, kij_tune=None):
    """ Model = SRK, PR76, PR78, RK """
    c = components
    temp = 350
    count = 0
    if model==None or model == "Ideal":
        while True:
            
            def f(temp):
                P_sum = 0
                for i in range(0,len(c)):
                    P_sum += z[i]*P / T_prop.general( c[i].VaporPressure, temp )
                return P_sum - 1
            
            def df():
                return (f(temp + 0.0001) - f(temp -0.0001))/0.0002
            
            f_temp = f(temp)
            df_temp = df()
            
            if abs(f_temp) < 0.0001:
                break
            else:
                temp = temp - f_temp/df_temp
            if count > 100:
                raise Warning("Dew point calculation not converged in 100 iterations!")
        return temp
    else:
        vapor_method = {"SRK":eos.SRK.phi_vapor,"Ideal":eos.Ideal.phi_vapor,"PR76":eos.PR76.phi_vapor,"PR78":eos.PR76.phi_vapor,"RK":eos.RK.phi_vapor} 
        while True:
            
            def f(temp):
                P_sum = 0
                phi = vapor_method[model](c,temp,P,z, kij_input, kij_tune)
                for i in range(0,len(c)):
                    P_sum += phi[i]*z[i]*P / T_prop.general( c[i].VaporPressure, temp )
                return P_sum - 1
            
            def df():
                return (f(temp + 0.0001) - f(temp -0.0001))/0.0002
            
            f_temp = f(temp)
            df_temp = df()
            
            if abs(f_temp) < 0.0001:
                break
            else:
                temp = temp - f_temp/df_temp
            if count > 100:
                raise Warning("Dew point calculation not converged in 100 iterations!")
        return temp

def PT_flash (components, pressure, temperature, fractions, models, mode=None, kij_input=None, kij_tune=None):
    """ Pressure-Temperature flash for VLE system.
    :param components: array containing chemical class objects
    :param models: 1:Vapor, 2:Liquid K value methods/models
    :param mode: 'ideal', 'gamma-phi', 'eos' 
    returns V/F, x, y"""
    
    cs = components #Chemical object array
    P = pressure
    T = temperature
    z = fractions
    check_fractions = round(sum(z),2)
    if check_fractions != 1.0:
        for i in range(0,len(z)):
            z[i] = z[i]/check_fractions
        print("WARNING! Sum of the entered molar component fractions isn't equal to 1. Molar fractions are normalized as: ",z)
    
    lower_limit = bubble_T(cs, z, P, models[1], kij_input, kij_tune)
    upper_limit = dew_T(cs, z, P, models[0], kij_input, kij_tune)
    
    if T < lower_limit or T > upper_limit:
        if T < lower_limit:
            return 0, z, [0,0]
        elif T > upper_limit:
            return 1, [0,0], z

    else:
        def beta(K):
            if len(K) == 2:# Binary mixture
                K1 = K[0]; K2 = K[1]; z1 = z[0]; z2 = z[1]
                guess = (-K1*z1 - K2*z2 + z1 + z2)/(K1*K2*z1 + K1*K2*z2 - K1*z1 - K1*z2 - K2*z1 - K2*z2 + z1 + z2)
            else:
                eps = 1e-3
                def f(guess,i):
                    f_beta = (K[i] - 1)*z[i] / (1 + guess*(K[i]-1))
                    return f_beta
                def df(guess,i):
                    df_beta = (-z[i]*(K[i]-1)**2) / (1+guess*(K[i]-1))**2
                    return df_beta
                guess = 0.5
                count = 0

                while True:
                    f_beta = 0; df_beta = 0; 
                    for i in range(0,len(cs)):
                        f_beta += f(guess,i)
                        df_beta += df(guess,i)
                    beta_new = guess - f_beta/df_beta
                    if abs((beta_new-guess)/guess) < eps:
                        break
                    guess = beta_new
                    count += 1
                    if count > 100:
                        raise Exception("Beta function not converged!")
            return guess

        #---Calculation of initial values---
        #Ki estimation using Wilson equation
        K = [] # K values holder
        x = [] # initial x fractions
        y = [] # initial y fractions
        
        for i in range(0,len(cs)):
            #Ki = Pc/P * exp(5.37*(1+w)*(1-Tc/T))
            Pc = float(cs[i].CriticalPressure)
            Tc = float(cs[i].CriticalTemperature)
            w = float(cs[i].AcentricityFactor)
            K.append((Pc/P)*exp(5.37*(1+w)*(1-Tc/T)))
        #K = [1.1, 1.1]
        beta_initial = beta(K)
        for t in range(0,len(cs)):    
            x.append( z[t]/(1+beta_initial*(K[t]-1)) )
            y.append( x[t]*K[t] )

        epsilon = 1e-7; counter = 0
        fug_V = []; fug_L = []
        V = beta_initial

        while True:
            vapor_method = {"SRK":eos.SRK.phi_vapor,"Ideal":eos.Ideal.phi_vapor,"PR76":eos.PR76.phi_vapor,"PR78":eos.PR76.phi_vapor,"RK":eos.RK.phi_vapor} 
            liquid_method = {"Ideal":activity.Ideal.gamma, "NRTL":activity.NRTL.gamma,"Uniquac":activity.Uniquac.gamma,"Unifac":activity.Unifac.gamma,"Dortmund":activity.Dortmund.gamma}  
            control = 0 # Control statement for convergence    
            
            if mode == "Ideal" or mode == None:
                for i in range(0,len(cs)):
                    Psat = T_prop.general(cs[i].VaporPressure, T)
                    fug_V.append( y[i]*P )
                    fug_L.append( x[i]*Psat )
                
            elif mode == "gamma-phi":
                phi = vapor_method[models[0]](cs,T,P,y, kij_input, kij_tune) 
                gamma = liquid_method[models[1]](cs,T,x)
                
                for i in range(0,len(cs)):
                    Psat = T_prop.general(cs[i].VaporPressure, T)
                    fug_V.append( phi[i]*y[i]*P )
                    fug_L.append( gamma[i]*x[i]*Psat )
            
            elif mode == "eos":
                phi_v = vapor_method[models[0]](cs,T,P,y, kij_input, kij_tune)  
                phi_l = liquid_method[models[0]](cs,T,P,x, kij_input, kij_tune) 
                for i in range(0,len(cs)):
                    fug_V.append( phi_v[i]*y[i])
                    fug_L.append( phi_l[i]*x[i])
                
            for k in range(0,len(cs)):
                control += y[k]*log(fug_L[k] / fug_V[k]) 
            
            counter+=1
            if abs(control) < epsilon:
                if V > 1.0001:
                    V = 1; y = z; x = [0,0]
                    
                elif V < 0.0001:
                    V = 0; x = z; y = [0,0]
                break
            
            else:
                for j in range(0,len(K)):
                    K[j] =  K[j]*fug_L[j]/fug_V[j] #fug_L[j]/fug_V[j]*y[j]/x[j]
                V = beta(K) #Update vapor fraction 
                
                for t in range(len(K)):
                    x[t] = z[t]/(1 + V*(K[t]-1))
                    y[t] = K[t]*x[t]
                fug_V = []; fug_L = []
            
            if counter > 100:
                raise Exception("Main loop in PT flash couldn't converged in 100 loops.")
                break
        return round(V,5), [round(i,5) for i in x], [round(i,5) for i in y]


"""
****EXAMPLE****

chem1 = get_chemical("1921")
chem2 = get_chemical("1102")

PT_flash([chem1,chem2],101325,352.5,[0.35, 0.65], models=["SRK","NRTL"],mode= "gamma-phi")

"""
