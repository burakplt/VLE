"""Activity coefficient methods"""
from math import log as ln
from math import exp, sqrt
import time

class Ideal():
    def gamma(components, temp, fractions):

        gammas = []
        for i in range(0, len(components)):
            gammas.append(1)
        return gammas

class Uniquac():
    """UNIQUAC model based activity coefficient calculations."""

    def gamma(components,temperature,fractions):
        
        cs = components 
        T = temperature
        x = fractions
        for item in x:
            if item == 0:
                item = 1E-05
              
        r = []; q = []; qp = []
        for k in range(0,len(cs)):
            r.append( float(cs[k].UniquacR) )
            q.append( float(cs[k].UniquacQ) )
            qp.append( float(cs[k].UniquacQP) )
        
        #---Calculate teta and fi values for each substance
        teta = []; fi = []; tetap = [] #tetap = teta' prime value for gR calc.
        for i in range(0,len(cs)):
            fi_nom = x[i]*r[i] #fi Nominator
            fi_denom = 0  #fi Denominator
            teta_nom = x[i]*q[i]
            teta_denom = 0
            tetap_nom = x[i]*qp[i]
            tetap_denom = 0
            for j in range(0,len(cs)):
                fi_denom += x[j]*r[j]
                teta_denom += x[j]*q[j]
                tetap_denom += x[j]*qp[j]
            fi.append(fi_nom/fi_denom) #Fi value of the i. component
            teta.append(teta_nom/teta_denom)  #teta value of the i. component
            tetap.append(tetap_nom/tetap_denom) #teta' prime value of the i. component
        
        def a_ij(id1, id2):
            file_path = "Models\\uniquac.txt"
            with open(file_path, 'r') as f:
                isFound = False # Is parameters found?
                for line in f.readlines():
                    aux = line.split(';')
                    if aux[0] == id1 and aux[1] == id2:
                        a12 = aux[2]
                        isFound = True
                    elif aux[0] == id2 and aux[1] == id1:
                        a12 = aux[3]
                        isFound = True
            if isFound:
                return float(a12) #units.mol_enthalpy(float(a12),"CGS","SI") #Convert to kJ/kmol
            else: 
                print('No parameters were found!')

        def tau(i,j):
            """Calculates tau_ij values"""
            if i == j:
                return 1
            else:
                id1 = cs[i].LibraryIndex
                id2 = cs[j].LibraryIndex
                return exp( -a_ij(id1,id2)/(1.9872*T)) #R = 1.9872 cal/mol.K 
        
        taus = {}
        for i in range(0,len(cs)):
            for j in range(0,len(cs)):
                taus[(i,j)] = tau(i,j)

        def unsymmetric():
            l = [0,0]
            l[0] = 5*(r[0]-q[0]) - (r[0]-1)
            l[1] = 5*(r[1]-q[1]) - (r[1]-1)
            C1 = ln(fi[0]/x[0]) + 5*q[0]*ln(teta[0]/fi[0]) + fi[1]*(l[0]-r[0]*l[1]/r[1]) 
            R1 = qp[0]*ln(tetap[0]+tetap[1]*taus[(1,0)]) + tetap[1]*qp[0]*(taus[(1,0)]/(tetap[0]+tetap[1]*taus[(1,0)]) - taus[(0,1)]/(tetap[1]+tetap[0]*taus[(0,1)]))  
            return exp(C1+R1)
        
        def symmetric():
            C = []; R = []
            for i in range(0,len(cs)):
                C.append( 1 + ln(fi[i]/x[i]) - fi[i]/x[i] -5*q[i]*( 1+ ln(fi[i]/teta[i])- fi[i]/teta[i] ))
                for j in range(0,len(cs)):
                    if i != j :
                        R.append( qp[i]*( 1- ln( tetap[j]*taus[(j,i)]+tetap[i] )- tetap[j]*taus[(i,j)]/(tetap[j]+tetap[i]*taus[(i,j)]) - tetap[i]/(tetap[j]*taus[(j,i)]+tetap[i]) ) )            
            return exp(C[0]+R[0]), exp(C[1]+R[1])
        
        return symmetric()
        
class NRTL():
    """NRTL Activity coefficient calculations"""

    def gamma(components, temperature, fractions):

        cs = components
        T = temperature
        x = fractions

        def a_ij(id1, id2):
            file_path = "Models\\nrtl.txt"
            with open(file_path, 'r') as f:
                isFound = False # Is parameters found?
                for line in f.readlines():
                    aux = line.split(';')
                    if aux[0] == id1 and aux[1] == id2:
                        a12 = aux[2]
                        alfa = aux[4]
                        isFound = True
                    elif aux[0] == id2 and aux[1] == id1:
                        alfa = aux[4]
                        a12 = aux[3]
                        isFound = True
            if isFound:
                return float(a12), float(alfa) #units.mol_enthalpy(float(a12),"CGS","SI") #Convert to kJ/kmol
            else: 
                print('WARNING!: No parameters were found for a_ij! Default parameters were used')
                return 100, 0.5 #Default parameters 

        aij = {}
        for i in range(0,len(cs)):
            for j in range(0,len(cs)):
                if i != j:
                    aij[(i,j)] = a_ij( cs[i].LibraryIndex, cs[j].LibraryIndex )
        
        def tau(i,j):
            """Calculates tau_ij values"""
            if i == j:
                return 0, 1
            else:
                aux_tau = aij[(i,j)]
                return aux_tau[0]/(1.9872*T), aux_tau[1]  #R = 1.9872 cal/mol.K
            
        def G(i,j):
            """Calculates Gij value"""
            if i == j:
                return 1
            else:
                aux_G = tau(i,j)
                return exp(-aux_G[1] * aux_G[0] )

        S = []; C = []
        for i in range (0, len(cs)):
            aux1 = 0; aux2 = 0
            for j in range(0,len(cs)):
                aux1 += x[j]*G(j,i) 
                aux2 += x[j]*G(j,i) * tau(j,i)[0]
            S.append(aux1)
            C.append(aux2)
        
        gamma = []
        for i in range(0,len(cs)):
            aux_k = 0
            for k in range(0,len(cs)):
                aux_k += x[k]*G(i,k)*(tau(i,k)[0] - C[k]/S[k])/S[k]
            gamma.append( exp( C[i]/S[i] + aux_k) )

        return gamma
    
    def gamma2(components, temperature, fractions):
        """Activity coefficients for binary mixture"""
        cs = components
        T = temperature
        x = fractions

        def a_ij(id1, id2):
            file_path = "Models\\nrtl.txt"
            with open(file_path, 'r') as f:
                isFound = False # Is parameters found?
                for line in f.readlines():
                    aux = line.split(';')
                    if aux[0] == id1 and aux[1] == id2:
                        a12 = aux[2]
                        alfa = aux[4]
                        isFound = True
                    elif aux[0] == id2 and aux[1] == id1:
                        alfa = aux[4]
                        a12 = aux[3]
                        isFound = True
            if isFound:
                return float(a12), float(alfa) #units.mol_enthalpy(float(a12),"CGS","SI") #Convert to kJ/kmol
            else: 
                print('WARNING!: No parameters were found for a_ij! Default parameters were used')
                return 100, 0.5 #Default parameters 

        aij = {}
        for i in range(0,len(cs)):
            for j in range(0,len(cs)):
                aij[(i,j)] = a_ij(cs[i].LibraryIndex, cs[j].LibraryIndex )

        def tau(i,j):
            """Calculates tau_ij values"""
            if i == j:
                return 0, 1
            else:
                aux_tau = aij[(i,j)]
                return aux_tau[0]/(1.9872*T), aux_tau[1]  #R = 1.9872 cal/mol.K
            
        def G(i,j):
            """Calculates Gij value"""
            if i == j:
                return 1
            else:
                aux_G = tau(i,j)
                return exp(-aux_G[1] * aux_G[0] )

        gama1 = (x[1]**2)*( tau(1,0)[0]*G(1,0)**2 /(x[0]+x[1]*G(1,0))**2 + tau(0,1)[0]*G(0,1)/ (x[1]+x[0]*G(0,1))**2 )
        gama2 = (x[0]**2)*( tau(0,1)[0]*G(0,1)**2 /(x[1]+x[0]*G(0,1))**2 + tau(1,0)[0]*G(1,0)/ (x[0]+x[1]*G(1,0))**2 )
        
        return exp(gama1), exp(gama2)

class Dortmund():
    """Modified Unifac Dortmund model"""
    def gamma(components,temperature,fractions):
        
        cs = components 
        T = temperature
        x = fractions
        for item in x:
            if item == 0:
                item = 1E-05
        
        # Get Q and R values for groups
        groupi = []; groupk = {}; ip = {}
        file_path = "Models\\modfac.txt"
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for i in range(0,len(cs)):
                groups = cs[i].ModifiedUnifac
                rk_data = []
                for pair in groups:
                    for line in lines:
                        aux = line.split(';')
                        if aux[3] == str(pair[0]):
                            ip[pair[0]] = int(aux[0]) 
                            if pair[0] in groupk.keys():
                                groupk[pair[0]][0].append((i,pair[1]))
                            else:
                                groupk[pair[0]] = ([(i, pair[1])], float(aux[4]), float(aux[5]))
                            rk_data.append( (pair[0], pair[1], float(aux[4]), float(aux[5])) )
                            break 
                groupi.append(rk_data)               
        #groupk= {17: ([(0, 1)], 0.92, 1.4), 1: ([(1, 1)], 0.9011, 0.848), 2: ([(1, 1)], 0.6744, 0.54), 15: ([(1, 1)], 1.0, 1.2)}
        
        #Calculate r and q values for components
        r = []; q = []
        for i in range(0,len(cs)):
            ri = 0; qi = 0
            for data in groupi[i]:
                ri += data[1]*data[2]
                qi += data[1]*data[3]
            r.append(ri)
            q.append(qi)        
        
        # Calculation of residual and combinatorial parts
        # ln gamma_k = Qk*[ 1-ln(sum(tetai*taui,k)) - sum [ (tetai*taui,m)/sum(tetaj*tauj,m)]
        # Calculate activity coefficients for each group
        
        group_names = [] # Get group numbers 
        for key in groupk.keys():
            group_names.append(key)
        
        def X(k):
            """Calculates group fraction for k"""
            aux_group = groupk[k] 
            aux1 = 0; aux2 = 0
            for item in aux_group[0]: #Item = (i, vi)
                vk = item[1]; i = item[0]
                aux1 += vk*x[i]
            
            for index in group_names:
                aux_grp = groupk[index][0]
                for itm in aux_grp:
                    aux2 += x[itm[0]]*itm[1]
            return aux1/aux2
        
        def tau(m,n):
            if m == n:
                return 1
            else:
                file_name = "Models\\modfac_ip.txt"
                found = False
                m = ip[m]; n = ip[n]
                with open(file_name, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        line = line.split()

                        if m == n:
                            aij = 0; bij = 0; cij = 0
                            found = True
                            break
                        elif int(line[0]) == m and int(line[1]) == n:
                            aij = float(line[2])
                            found = True
                            bij = float(line[3])
                            cij = float(line[4])
                            break
                        elif int(line[0]) == n and int(line[1]) == m:
                            aij = float(line[5])
                            found = True
                            bij = float(line[6])
                            cij = float(line[7])
                            break
                if found:
                    return exp( -(aij + bij*T + cij*T**2)/T)
                else:
                    print("WARNING! No MODFAC interaction parameters were found for groups",m,n)
                    return exp(-50/T) #default value
        
        taus = {}
        for m in group_names:
            for n in group_names:
                taus[(m,n)] = tau(m,n)

        Xk = [] #Calculate and store Xk values
        for k in group_names:
            Xk.append(X(k))
        
        Xi = [] #Calculate and store Xk values for pure components
        def X2(k, xi):
            """Calculates group fraction for k"""
            aux_group = groupk[k] 
            aux1 = 0; aux2 = 0
            for item in aux_group[0]: #Item = (i, vi)
                vk = item[1]; i = item[0]
                aux1 += vk*xi[i]
            for index in group_names:
                aux_grp = groupk[index][0]
                for itm in aux_grp:
                    aux2 += xi[itm[0]]*itm[1]
            return aux1/aux2

        def teta(k):
            """Teta value for group m"""
            Qk = groupk[k][2]
            kk = group_names.index(k)
            aux = 0
            for n in group_names:
                nk = group_names.index(n)
                Qn = groupk[n][2]
                aux += Qn*Xk[nk]
            tet = groupk[k][2]*Xk[kk]/aux
            return tet

        for i in range(0,len(cs)):
            ki = []
            for k in group_names:# TODO loop sırasını değiştir i dışa k içe
                xi = x.copy()
                for j in range(0,len(xi)): 
                    if i==j:
                        xi[j] = 1
                    else:
                        xi[j] = 0
                ki.append( X2(k,xi) ) 
            Xi.append(ki)
        
        def tetai(k, i):
            """Teta value for group m in pure component"""
            Qk = groupk[k][2]
            kk = group_names.index(k)
            aux = 0
            for n in group_names:
                nk = group_names.index(n)
                Qn = groupk[n][2]
                aux += Qn*Xi[i][nk]
            teti = groupk[k][2]*Xi[i][kk]/aux
            return teti
        
        teta_k = []; teta_ki = []
        for i in range(0,len(cs)):
            pure_k = []
            for k in group_names:
                pure_k.append(tetai(k,i))
            teta_ki.append(pure_k)
        
        for k in group_names:
                teta_k.append(teta(k))
        
        activity_R = [] #Residual part for activity coefficient ln gammaR
        for i in range(0,len(cs)):
            ln_gamma_R = 0
            
            for k in group_names:
                vk = 0
                for t in groupk[k][0]:
                    if t[0] == i:
                        vk = t[1]
                
                Qk = groupk[k][2]
                kk = group_names.index(k)
                nom = 0; aux = 0
                nom_i = 0; aux_i = 0
                
                for m in group_names:
                    denom_i = 0; denom = 0
                    mm = group_names.index(m)
                    for n in group_names:
                        nn = group_names.index(n)
                        denom += teta_k[nn]*taus[(n,m)]
                        denom_i += teta_ki[i][nn]*taus[(n,m)]
                        
                    nom += teta_k[mm]*taus[(k,m)]/denom
                    aux += teta_k[mm]*taus[(m,k)]
                    nom_i += teta_ki[i][mm]*taus[(k,m)]/denom_i
                    aux_i += teta_ki[i][mm]*taus[(m,k)]
                
                ln_gamma_k = Qk*(1- ln(aux) - nom )
                ln_gamma_ki = Qk*(1- ln(aux_i) - nom_i )
                ln_gamma_R += vk*(ln_gamma_k - ln_gamma_ki)
                #print(k, ln_gamma_k)
            activity_R.append(ln_gamma_R)
        
        activity_C = []
        #Gamma combinatorial for components
        V = []; F = []; Vp = [] # V' modified dortmund
        for i in range (0,len(cs)):
            aux_r = 0; aux_q = 0; aux_rp = 0
            for j in range(0,len(cs)):
                aux_r += r[j]*x[j]
                aux_rp += (r[j]**0.75)*x[j]
                aux_q += q[j]*x[j]
            V.append(r[i]/aux_r)
            Vp.append((r[i]**0.75)/aux_rp)
            F.append(q[i]/aux_q)
        
        for i in range(0, len(cs)):
            aux = 1 - Vp[i]+ ln(Vp[i]) - 5*q[i]*( 1- V[i]/F[i]+ ln(V[i]/F[i]) )
            activity_C.append(aux)
        
        activity_coefficients = []
        for i in range(0,len(cs)):
            activity_coefficients.append( exp(activity_C[i] + activity_R[i]) )
        
        return activity_coefficients

class Unifac():
    """Unifac model activity coefficient"""
    def gamma(components,temperature,fractions):
        
        cs = components 
        T = temperature
        x = fractions
        for item in x:
            if item == 0:
                item = 1E-05
        
        # Get Q and R values for groups
        groupi = []; groupk = {}; ip = {}
        file_path = "Models\\unifac.txt"
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for i in range(0,len(cs)):
                groups = cs[i].UnifacVLE
                rk_data = []
                for pair in groups:
                    for line in lines:
                        aux = line.split(',')
                        if aux[1] == str(pair[0]):
                            ip[pair[0]] = int(aux[0]) 
                            if pair[0] in groupk.keys():
                                groupk[pair[0]][0].append((i,pair[1]))
                            else:
                                groupk[pair[0]] = ([(i, pair[1])], float(aux[4]), float(aux[5]))
                            rk_data.append( (pair[0], pair[1], float(aux[4]), float(aux[5])) )
                            break 
                groupi.append(rk_data)               
        #groupk= {17: ([(0, 1)], 0.92, 1.4), 1: ([(1, 1)], 0.9011, 0.848), 2: ([(1, 1)], 0.6744, 0.54), 15: ([(1, 1)], 1.0, 1.2)}
        
        #Calculate r and q values for components
        r = []; q = []
        for i in range(0,len(cs)):
            ri = 0; qi = 0
            for data in groupi[i]:
                ri += data[1]*data[2]
                qi += data[1]*data[3]
            r.append(ri)
            q.append(qi)        
       
        # Calculation of residual and combinatorial parts
        # ln gamma_k = Qk*[ 1-ln(sum(tetai*taui,k)) - sum [ (tetai*taui,m)/sum(tetaj*tauj,m)]
        # Calculate activity coefficients for each group
        
        group_names = [] # Get group numbers 
        for key in groupk.keys():
            group_names.append(key)
        
        def X(k):
            """Calculates group fraction for k"""
            aux_group = groupk[k] 
            aux1 = 0; aux2 = 0
            for item in aux_group[0]: #Item = (i, vi)
                vk = item[1]; i = item[0]
                aux1 += vk*x[i]
            
            for index in group_names:
                aux_grp = groupk[index][0]
                for itm in aux_grp:
                    aux2 += x[itm[0]]*itm[1]
            return aux1/aux2
        
        def tau(m,n):
            if m == n:
                return 1
            else:
                file_name = "Models\\unifac_ip.txt"
                found = False
                m = ip[m]; n = ip[n]
                with open(file_name, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        line = line.split("\t")
                        if int(line[0]) == m and int(line[2]) == n:
                            aij = float(line[4]); found = True
                        elif int(line[0]) == m and int(line[2]) == n:
                            aij = float(line[5]); found = True
                if found:
                    return exp(-aij/T)
                else:
                    print("WARNING! No UNIFAC interaction parameters were found for groups",m,n)
                    return exp(-50/T) #default value
        
        taus = {}
        for m in group_names:
            for n in group_names:
                taus[(m,n)] = tau(m,n)

        Xk = [] #Calculate and store Xk values
        for k in group_names:
            Xk.append(X(k))
        
        Xi = [] #Calculate and store Xk values for pure components
        def X2(k, xi):
            """Calculates group fraction for k"""
            aux_group = groupk[k] 
            aux1 = 0; aux2 = 0
            for item in aux_group[0]: #Item = (i, vi)
                vk = item[1]; i = item[0]
                aux1 += vk*xi[i]
            for index in group_names:
                aux_grp = groupk[index][0]
                for itm in aux_grp:
                    aux2 += xi[itm[0]]*itm[1]
            return aux1/aux2

        def teta(k):
            """Teta value for group m"""
            Qk = groupk[k][2]
            kk = group_names.index(k)
            aux = 0
            for n in group_names:
                nk = group_names.index(n)
                Qn = groupk[n][2]
                aux += Qn*Xk[nk]
            tet = groupk[k][2]*Xk[kk]/aux
            return tet

        t5 = time.process_time()
        for i in range(0,len(cs)):
            ki = []
            for k in group_names:
                xi = x.copy()
                for j in range(0,len(xi)): 
                    if i==j:
                        xi[j] = 1
                    else:
                        xi[j] = 0
                ki.append( X2(k,xi) ) 
            Xi.append(ki)

        def tetai(k, i):
            """Teta value for group m in pure component"""
            Qk = groupk[k][2]
            kk = group_names.index(k)
            aux = 0
            for n in group_names:
                nk = group_names.index(n)
                Qn = groupk[n][2]
                aux += Qn*Xi[i][nk]
            teti = groupk[k][2]*Xi[i][kk]/aux
            return teti
        
        teta_k = []; teta_ki = []
        for i in range(0,len(cs)):
            pure_k = []
            for k in group_names:
                pure_k.append(tetai(k,i))
            teta_ki.append(pure_k)
        
        for k in group_names:
                teta_k.append(teta(k))
        
        activity_R = [] #Residual part for activity coefficient ln gammaR
        for i in range(0,len(cs)):
            ln_gamma_R = 0
            
            for k in group_names:
                vk = 0
                for t in groupk[k][0]:
                    if t[0] == i:
                        vk = t[1]
                
                Qk = groupk[k][2]
                kk = group_names.index(k)
                nom = 0; aux = 0
                nom_i = 0; aux_i = 0
                
                for m in group_names:
                    denom_i = 0; denom = 0
                    mm = group_names.index(m)
                    for n in group_names:
                        nn = group_names.index(n)
                        denom += teta_k[nn]*taus[(n,m)]
                        denom_i += teta_ki[i][nn]*taus[(n,m)]
                        
                    nom += teta_k[mm]*taus[(k,m)]/denom
                    aux += teta_k[mm]*taus[(m,k)]
                    nom_i += teta_ki[i][mm]*taus[(k,m)]/denom_i
                    aux_i += teta_ki[i][mm]*taus[(m,k)]
                
                ln_gamma_k = Qk*(1- ln(aux) - nom )
                ln_gamma_ki = Qk*(1- ln(aux_i) - nom_i )
                ln_gamma_R += vk*(ln_gamma_k - ln_gamma_ki)
                
            activity_R.append(ln_gamma_R)
        
        activity_C = []
        #Gamma combinatorial for components
        V = []; F = []
        for i in range (0,len(cs)):
            aux_r = 0; aux_q = 0
            for j in range(0,len(cs)):
                aux_r += r[j]*x[j]
                aux_q += q[j]*x[j]
            V.append(r[i]/aux_r)
            F.append(q[i]/aux_q)
        
        for i in range(0, len(cs)):
            aux = 1 - V[i]+ ln(V[i]) - 5*q[i]*( 1- V[i]/F[i]+ ln(V[i]/F[i]) )
            activity_C.append(aux)

        activity_coefficients = []
        for i in range(0,len(cs)):
            activity_coefficients.append( exp(activity_C[i] + activity_R[i]) )
        
        return activity_coefficients

