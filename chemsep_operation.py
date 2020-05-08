"""ChemSep database operations. Gather serialized XML data,
    convert paramaters and add functions"""
import math
import numpy as np 
import xml.etree.ElementTree as ET
from pickle import load, dump


class Chemical(object):
    
    def __init__(self, name, lib_index):
        self.name = name
        self.LibraryIndex = lib_index
        """
        self.CriticalTemperature = "K"
        self.CriticalPressure = "Pa"
        self.CriticalVolume = "m3"
        self.CriticalCompressibility = ""
        self.NormalBoilingPointTemperature = ""
        self.NormalMeltingPointTemperature = ""
        self.TriplePointTemperature = ""
        self.TriplePointPressure = ""
        self.MolecularWeight = ""
        self.LiquidVolumeAtNormalBoilingPoint = ""
        self.AcentricityFactor = ""
        self.SolubilityParameter = ""
        self.DipoleMoment = ""
        self.HeatOfFormation = ""
        self.GibbsEnergyOfFormation = ""
        self.AbsEntropy = ""
        self.HeatOfFusionAtMeltingPoint = ""
        self.HeatOfCombustion = ""
        self.COSTALDVolume = ""
        self.DiameterLJ = ""
        self.EnergyLJ = ""
        self.RacketParameter = ""
        self.FullerVolume = ""
        self.Parachor = ""
        self.SpecificGravity = ""
        self.CostaldAcentricFactor = ""
        self.WilsonVolume = ""
        self.ChaoSeaderAcentricFactor = ""
        self.ChaoSeaderSolubilityParameter = ""
        self.ChaoSeaderLiquidVolume = ""
        self.MatthiasCopemanC1 = ""
        self.MatthiasCopemanC2 = ""
        self.MatthiasCopemanC3 = ""
        self.UniquacR = ""
        self.UniquacQ = ""
        self.UniquacQP = ""
        self.ApiSrkS1 = ""
        self.ApiSrkS2 = ""
        self.UnifacVLE = [] #Array
        self.UnifacLLE = [] #Array
        self.Asog = [] #Array
        self.GCmethod = [] #Array
        self.Umr = [] #Array
        self.ModifiedUnifac = [] #Array """


def eq1 (A):
    """Chemsep equation 1
    :param A: Equation parameter"""
    return A

def eq2 (A, B, T):
    """Chemsep equation 2
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param T: Temperature in K"""
    return A + B*T

def eq3 (A, B, C, T):
    """Chemsep equation 3
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param T: Temperature in K"""
    return A + B*T + C*(T**2)

def eq4 (A, B, C, D, T):
    """Chemsep equation 4
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param T: Temperature in K"""
    return A + B*T + C*(T**2) + D*(T**3)

def eq5 (A, B, C, D, E, T):
    """Chemsep equation 5
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param E: Equation parameter E
    :param T: Temperature in K"""
    return A + B*T + C*(T**2) + D*(T**3) + E*(T**4)

def eq10 (A, B, C, T):
    """Chemsep equation 10
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param T: Temperature in K"""
    return math.exp( A - B/(T+C) )

def eq11 (A):
    """Chemsep equation 11
    :param A: Equation parameter A"""
    return math.exp(A)

def eq12 (A, B, T):
    """Chemsep equation 12
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param T: Temperature in K"""
    return math.exp(A + B*T)

def eq13 (A, B, C, T):
    """Chemsep equation 13
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param T: Temperature in K"""
    return math.exp( A + B*T + C*(T**2) )

def eq14 (A, B, C, D, T):
    """Chemsep equation 14
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param T: Temperature in K"""
    return math.exp( A + B*T + C*(T**2) + D*(T**3) )

def eq15 (A, B, C, D, E, T):
    """Chemsep equation 15
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param E: Equation parameter E
    :param T: Temperature in K"""
    return math.exp( A + B*T + C*(T**2) + D*(T**3) + E*(T**4) )

def eq16 (A, B, C, D, E, T):
    """Chemsep equation 16
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param E: Equation parameter E
    :param T: Temperature in K"""
    return A + math.exp( B/T + C + D*T + E*(T**2) )

def eq17 (A, B, C, D, E, T):
    """Chemsep equation 17
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param E: Equation parameter E
    :param T: Temperature in K"""
    return A + math.exp( B + C*T + D*(T**2) + E*(T**3) )

def eq100 (A, B, C, D, E, T):
    """Chemsep equation 100
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param E: Equation parameter E
    :param T: Temperature in K"""
    return A + B*T + C*(T**2) + D*(T**3) + E*(T**4)

def eq101 (A, B, C, D, E, T):
    """Chemsep equation 101
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param E: Equation parameter E
    :param T: Temperature in K"""
    return math.exp( A + B/T + C*math.log(T) + D*(T**E) )

def eq102 (A, B, C, D, T):
    """Chemsep equation 102
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param T: Temperature in K"""
    return A*(T**B) / (1 + C/T + D/(T**2))

def eq104 (A, B, C, D, E, T):
    """Chemsep equation 104
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param E: Equation parameter E
    :param T: Temperature in K"""
    return A + B/T + C/(T**3) + D/(T**8) + E/(T**9)

def eq105 (A, B, C, D, T):
    """Chemsep equation 105
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param T: Temperature in K"""
    body = 1+(1-T/C)**D
    return A/math.pow(B,body)

def eq106 (A, B, C, D, E, Tc, T):
    """Chemsep equation 106
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param E: Equation parameter E
    :param Tc: Critical temperature
    :param T: Temperature in K"""
    Tr = T/Tc #Critical temperature
    body = B + C*Tr + D*(Tr**2) + E*(Tr**3)
    return A * math.pow(1-Tr, body)

def eq120 (A, B, C, T):
    """Chemsep equation 120
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param T: Temperature in K"""
    return A - B/(T+C)

def eq121 (A, B, C, D, E, T):
    """Chemsep equation 121
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param D: Equation parameter D
    :param E: Equation parameter E
    :param T: Temperature in K"""
    return A + B/T + C*math.log(T) + D*(T**E)

def eq208 (A, B, C, T):
    """Chemsep equation 208, Antoine equation
    :param A: Equation parameter A
    :param B: Equation parameter B
    :param C: Equation parameter C
    :param T: Temperature in K"""
    body = A - B/(T+C)
    return math.pow(10, body)

#Temperature correlated data graphs
#param = equation_parameters(chemsep_xml,2,2,50)
#equation_to_array(param[0],param[1],param[2],param[3])
def plot_property (eqno, ID, p, data_points):
    """Converts temperature dependent equations
    to x and y values. Plot the graph if needed.
    :param eqno:[int] Which equation is used?
    :param id: id[0]=name, id[1]=CAS-No
    :param p: Dict object with all parameters
    :p[0]=A, [1]=B [2,3...]=C,D... 
    :param data_points: Tmin, Tmax, data point number
    :return A list containing x,y tuples."""
    y = [] # f(x) values
    Tmin = data_points[0] # Lower temperature limit
    Tmax = data_points[1] # Higher temperature limit
    data = data_points[2] # Number of points

    if eqno == 1:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq1(params[0]) )
    
    elif eqno == 2:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq2(p[0], p[1], x[i]) )

    elif eqno == 3:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq3(p[0],p[1],p[2],x[i]) )
    
    elif eqno == 4:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq4(p[0],p[1],p[2],p[3],x[i]) )
    
    elif eqno == 5:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq5(p[0],p[1],p[2],p[3],p[4],x[i]) )
    
    elif eqno == 10:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq10(p[0],p[1],p[2],x[i]) )
    
    elif eqno == 11:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq11(p[0]) )
    
    elif eqno == 12:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq12(p[0],p[1],x[i]) )
    
    elif eqno == 13:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq13(p[0],p[1],p[2],x[i]) )
    
    elif eqno == 14:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq14(p[0],p[1],p[2],p[3],x[i]) )
    
    elif eqno == 15:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq15(p[0],p[1],p[2],p[3],p[4],x[i]) )
    
    elif eqno == 16:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq16(p[0],p[1],p[2],p[3],p[4],x[i]) )
    
    elif eqno == 17:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq17(p[0],p[1],p[2],p[3],p[4],x[i]) )
    
    elif eqno == 100:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq100(p[0],p[1],p[2],p[3],p[4],x[i]) )
    
    elif eqno == 101:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq101(p[0],p[1],p[2],p[3],p[4],x[i]) )
    
    elif eqno == 102:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq102(p[0],p[1],p[2],p[3],x[i]) )
    
    elif eqno == 104:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq104(p[0],p[1],p[2],p[3],p[4],x[i]) )
    
    elif eqno == 105:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq105(p[0],p[1],p[2],p[3],x[i]) )
    
    elif eqno == 106:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq106(p[0],p[1],p[2],p[3],p[4],p[5],x[i]) )
    
    elif eqno == 120:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq120(p[0],p[1],p[2],x[i]) )
    
    elif eqno == 121:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq121(p[0],p[1],p[2],p[3],p[4],x[i]) )

    elif eqno == 208:
        x = np.linspace(Tmin, Tmax, data)
        for i in range(0, data):
            y.append( eq208(p[0],p[1],p[2],x[i]) )
    else:
        print("Invalid input data! Check the values.\n")
  
    
    plt.plot(x,y)
    plt.title(ID[1]+" - "+ID[2])
    plt.xlabel("Temperature (K)")
    plt.ylabel(ID[3])
    plt.show()

#Importer function from XML file all required parameters
def equation_parameters(file, index, prop_no, resolution):
    """Creates Temperature dependent property plots
    :param file: File name path
    :param index: Index for the root of tree
    :param prop_no: Which property? Select integer value."1":"LiquidDensity", "2":"VaporPressure","3":"HeatOfVaporization","4":"LiquidHeatCapacityCp", "5":"IdealGasHeatCapacityCp", "6":"SecondVirialCoefficient","7":"LiquidViscosity", "8":"VaporViscosity", "9":"LiquidThermalConductivity","10":"VaporThermalConductivity", "11":"RPPHeatCapacityCp","12":"RelativeStaticPermittivity","13":"AntoineVaporPressure", "14":"LiquidViscosityRPS","15":"SolidDensity","16":"SurfaceTension","17":"SolidHeatCapacityCp"
    :param resplution: How many data points will be used?
    :return eqno, ID, p, data_points for the equation to array()"""

    tree = ET.parse(file)
    root = tree.getroot()
    properties = {"1":"LiquidDensity", "2":"VaporPressure","3":"HeatOfVaporization",
    "4":"LiquidHeatCapacityCp", "5":"IdealGasHeatCapacityCp", "6":"SecondVirialCoefficient",
    "7":"LiquidViscosity", "8":"VaporViscosity", "9":"LiquidThermalConductivity",
    "10":"VaporThermalConductivity", "11":"RPPHeatCapacityCp","12":"RelativeStaticPermittivity",
    "13":"AntoineVaporPressure", "14":"LiquidViscosityRPS","15":"SolidDensity","16":"SurfaceTension",
    "17":"SolidHeatCapacityCp"}

    child = root[index] # Compound 
    ID = [0,0,0,0] # Compound identity
    data_points = [0,0,resolution] #Tmin, Tmax, data point
    p = [] #Parameters A,B...
    ID[0] = child.find('CAS').get('value') #CAS-No of the chemical
    ID[1] = child.find('CompoundID').get('value') #Name of the chemical
    prop = properties.get(str(prop_no))

    compound = child.find(prop)
    ID[2] = compound.get('name')
    ID[3] = compound.get('units')
    for item in compound:
        if item.tag == "eqno":
            eqno = int(item.get('value'))
        elif item.tag =="Tmin":
            data_points[0] = float(item.get('value'))
        elif item.tag =="Tmax":
            data_points[1] = float(item.get('value'))
        else:
            p.append(float(item.get('value')))
            
    return eqno, ID, p, data_points


class EosInterface(object):

    def ig_enthalpy(data, temperature):
        eqno = int(data[0]) #Equation no 
        p = []
        a = 298.15
        b = temperature
        n = 9
        h = (b-a)/n
        for i in range(1,len(data)):
            p.append(float(data[i]))
        
        if eqno == 16:
            total = eq16(p[0],p[1],p[2],p[3],p[4], a) + eq16(p[0],p[1],p[2],p[3],p[4],b) #first and last terms
            for i in range(1,n):
                if (i%3 == 0):
                    total += 2 *eq16(p[0],p[1],p[2],p[3],p[4], a+i*h)
                else:
                    total += 3 * eq16(p[0],p[1],p[2],p[3],p[4], a+i*h)
            return 3*h*total/8
        
        elif eqno == 1:
            return p[0]
        
        elif eqno == 100:
            total = eq100(p[0],p[1],p[2],p[3],p[4], a) + eq100(p[0],p[1],p[2],p[3],p[4],b) #first and last terms
            for i in range(1,n):
                if (i%3 == 0):
                    total += 2 *eq100(p[0],p[1],p[2],p[3],p[4], a+i*h)
                else:
                    total += 3 * eq100(p[0],p[1],p[2],p[3],p[4], a+i*h)
            return 3*h*total/8
    
    def ig_entropy(data, temperature):
        eqno = int(data[0]) #Equation no 
        p = []
        a = 298.15
        b = temperature
        n = 9
        h = (b-a)/n
        for i in range(1,len(data)):
            p.append(float(data[i]))
        
        if eqno == 16:
            def eq16_S (A, B, C, D, E, T):
                return (A + math.exp( B/T + C + D*T + E*(T**2) ))/ T
            total = eq16_S(p[0],p[1],p[2],p[3],p[4], a) + eq16_S(p[0],p[1],p[2],p[3],p[4],b) #first and last terms
            for i in range(1,n):
                if (i%3 == 0):
                    total += 2 * eq16_S(p[0],p[1],p[2],p[3],p[4], a+i*h)
                else:
                    total += 3 * eq16_S(p[0],p[1],p[2],p[3],p[4], a+i*h)
            return 3*h*total/8
        
        elif eqno == 1:
            return p[0]
        
        elif eqno == 100:
            total = eq100(p[0],p[1],p[2],p[3],p[4], a) + eq100(p[0],p[1],p[2],p[3],p[4],b) #first and last terms
            for i in range(1,n):
                if (i%3 == 0):
                    total += 2 *eq100(p[0],p[1],p[2],p[3],p[4], a+i*h)
                else:
                    total += 3 * eq100(p[0],p[1],p[2],p[3],p[4], a+i*h)
            return 3*h*total/8
       
    def dh_vaporization(data, temperature, Tc):
        #eqno = int(data[0]) #Equation no 
        p = []
        a = 298.15 / float(Tc)
        b = temperature / float(Tc)
        n = 9
        h = (b-a)/n
        
        def eq_106(A,B,C,D,E,Tr):
            body = B + C*Tr + D*(Tr**2) + E*(Tr**3)
            return A * math.pow(1-Tr, body)
        
        for i in range(1,len(data)):
            p.append(float(data[i]))
        
        total = eq_106(p[0],p[1],p[2],p[3],p[4], a) + eq_106(p[0],p[1],p[2],p[3],p[4],b) #first and last terms
        for i in range(1,n):
            if (i%3 == 0):
                total += 2 * eq_106(p[0],p[1],p[2],p[3],p[4], a+i*h)
            else:
                total += 3 * eq_106(p[0],p[1],p[2],p[3],p[4], a+i*h)
        return 3*h*total/8

    def general(data, temperature):
        eqno = int(data[0])
        p = []
        for i in range(1,len(data)):
            p.append(float(data[i]))

        if eqno == 1:
            return eq1(p[0])
        
        elif eqno == 2:
            return eq2(p[0], p[1], temperature)

        elif eqno == 3:
            return eq3(p[0],p[1],p[2],temperature) 
        
        elif eqno == 4:
            return eq4(p[0],p[1],p[2],p[3],temperature)
        
        elif eqno == 5:
            return eq5(p[0],p[1],p[2],p[3],p[4],temperature)
        
        elif eqno == 10:
            return eq10(p[0],p[1],p[2],temperature)
        
        elif eqno == 11:
            return eq11(p[0])
        
        elif eqno == 12:
            return eq12(p[0],p[1],temperature)
        
        elif eqno == 13:
            return eq13(p[0],p[1],p[2],temperature)
        
        elif eqno == 14:
            return eq14(p[0],p[1],p[2],p[3],temperature)
        
        elif eqno == 15:
            return eq15(p[0],p[1],p[2],p[3],p[4],temperature)
        
        elif eqno == 16:
            return eq16(p[0],p[1],p[2],p[3],p[4],temperature)
        
        elif eqno == 17:
            return eq17(p[0],p[1],p[2],p[3],p[4],temperature)
        
        elif eqno == 100:
            return eq100(p[0],p[1],p[2],p[3],p[4],temperature)
        
        elif eqno == 101:
            return eq101(p[0],p[1],p[2],p[3],p[4],temperature)
        
        elif eqno == 102:
            return eq102(p[0],p[1],p[2],p[3],temperature)
        
        elif eqno == 104:
            return eq104(p[0],p[1],p[2],p[3],p[4],temperature)
        
        elif eqno == 105:
            return eq105(p[0],p[1],p[2],p[3],temperature)
        
        elif eqno == 106:
            return eq106(p[0],p[1],p[2],p[3],p[4],p[5],temperature)
        
        elif eqno == 120:
            return eq120(p[0],p[1],p[2],temperature)
        
        elif eqno == 121:
            return eq121(p[0],p[1],p[2],p[3],p[4],temperature)

        elif eqno == 208:
            return eq208(p[0],p[1],p[2],temperature)
        else:
            print("Invalid input data! Check the values.\n")

def create_chemical():
    """Search and serialize XML data into Chemical class object.
    """
    
    properties = ("CriticalTemperature","CriticalPressure","CriticalVolume","CriticalCompressibility",
    "NormalBoilingPointTemperature","NormalMeltingPointTemperature","TriplePointTemperature",
    "TriplePointPressure","MolecularWeight","LiquidVolumeAtNormalBoilingPoint","AcentricityFactor",
    "SolubilityParameter","DipoleMoment","HeatOfFormation","GibbsEnergyOfFormation","AbsEntropy",
    "HeatOfFusionAtMeltingPoint","HeatOfCombustion","COSTALDVolume","DiameterLJ",
    "EnergyLJ","RacketParameter","FullerVolume","Parachor","SpecificGravity",
    "CostaldAcentricFactor","WilsonVolume","ChaoSeaderAcentricFactor","ChaoSeaderSolubilityParameter",
    "ChaoSeaderLiquidVolume","MatthiasCopemanC1","MatthiasCopemanC2","MatthiasCopemanC3",
    "UniquacR","UniquacQ","UniquacQP","ApiSrkS1","ApiSrkS2","UnifacVLE","UnifacLLE",
    "Asog","GCmethod","Umr","ModifiedUnifac","CAS")
    
    T_properties = {"1":"LiquidDensity", "2":"VaporPressure","3":"HeatOfVaporization",
    "4":"LiquidHeatCapacityCp", "5":"IdealGasHeatCapacityCp", "6":"SecondVirialCoefficient",
    "7":"LiquidViscosity", "8":"VaporViscosity", "9":"LiquidThermalConductivity",
    "10":"VaporThermalConductivity", "11":"RPPHeatCapacityCp","12":"RelativeStaticPermittivity",
    "13":"AntoineVaporPressure", "14":"LiquidViscosityRPS","15":"SolidDensity","16":"SurfaceTension",
    "17":"SolidHeatCapacityCp"}
    
    chemsep_xml = 'Databases\\chemsep.xml'
    tree = ET.parse(chemsep_xml)
    root = tree.getroot()
    compound="" #Place holder for compound 

    for elem in root.iter(tag='compound'): #Find compound in XML file and locate
        
        compound  = elem
        name = compound[1].get('value')
        lib = compound[0].get('value')
        substance = Chemical(name,lib) #create a Chemical class object
        
        for item in properties:
            prop = compound.find(item, None)
            if prop != None:
                aux = [] # Placeholder for prop values
                for sub_prop in prop:
                    if sub_prop.tag == 'group':
                        aux.append( ( int(sub_prop.attrib['id']), int(sub_prop.attrib['value']) ) )
                if aux:
                    setattr(substance, prop.tag, aux)
                else:
                    setattr(substance, prop.tag, prop.get('value'))
        
        for key in T_properties.keys():
            T_prop = compound.find(T_properties.get(key), None)
            d = [] #Placeholder for data points
            if T_prop != None:
                for item in T_prop:
                    if item.tag == "eqno":
                        d.append(item.get('value'))
                    elif item.tag =="Tmin":
                        pass
                    elif item.tag =="Tmax":
                        pass
                    else:
                        d.append(item.get('value'))
                setattr(substance, T_prop.tag, d)
        
        file_name = lib
        pickle_out = open("chemicals\\"+file_name+".pickle", "wb")
        dump(substance, pickle_out)
        pickle_out.close()


def get_chemical(lib_index):
    
    path = "chemicals\\"+lib_index+".pickle"
    
    with open(path, "rb") as f:
        chemical = load(f)
    return chemical

