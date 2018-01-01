# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
You can freely use and distribute it. If you hava  

any problem, you could contact with us timely!

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com

Z: atomic number
L: principal quantum number
Zv: number of valence electrons
Rv: van der Waals atomic radius
Rc: covalent radius
m: atomic mass
V: van der Waals vloume
En: Sanderson electronegativity
alapha: atomic polarizability (10e-24 cm3)
IP: ionization potential (eV)
EA: electron affinity (eV)
"""

################################################################################
AtomProperty={
'H':{'Z':1,'L':1,'Zv':1,'Rv':1.17,'Rc':0.37,'m':1.01,'V':6.71,'En':2.59,'alapha':0.67,'IP':13.598,'EA':0.754},
'Li':{'Z':3,'L':2,'Zv':1,'Rv':1.82,'Rc':1.34,'m':6.94,'V':25.25,'En':0.89,'alapha':24.3,'IP':5.392,'EA':0.618},
'Be':{'Z':4,'L':2,'Zv':2,'Rv':0.0,'Rc':0.90,'m':9.01,'V':0.0,'En':1.81,'alapha':5.60,'IP':9.323,'EA':0.0},
'B':{'Z':5,'L':2,'Zv':3,'Rv':1.62,'Rc':0.82,'m':10.81,'V':17.88,'En':2.28,'alapha':3.03,'IP':8.298,'EA':0.277},
'C':{'Z':6,'L':2,'Zv':4,'Rv':1.75,'Rc':0.77,'m':12.01,'V':22.45,'En':2.75,'alapha':1.76,'IP':11.260,'EA':1.263},
'N':{'Z':7,'L':2,'Zv':5,'Rv':1.55,'Rc':0.75,'m':14.01,'V':15.60,'En':3.19,'alapha':1.10,'IP':14.534,'EA':0.0},
'O':{'Z':8,'L':2,'Zv':6,'Rv':1.40,'Rc':0.73,'m':16.00,'V':11.49,'En':3.65,'alapha':0.80,'IP':13.618,'EA':1.461},
'F':{'Z':9,'L':2,'Zv':7,'Rv':1.30,'Rc':0.71,'m':19.00,'V':9.20,'En':4.00,'alapha':0.56,'IP':17.423,'EA':3.401},
'Na':{'Z':11,'L':3,'Zv':1,'Rv':2.27,'Rc':1.54,'m':22.99,'V':49.00,'En':0.56,'alapha':23.6,'IP':5.139,'EA':0.548},
'Mg':{'Z':12,'L':3,'Zv':2,'Rv':1.73,'Rc':1.30,'m':24.31,'V':21.69,'En':1.32,'alapha':10.6,'IP':7.646,'EA':0.0},
'Al':{'Z':13,'L':3,'Zv':3,'Rv':2.06,'Rc':1.18,'m':26.98,'V':36.51,'En':1.71,'alapha':6.80,'IP':5.986,'EA':0.441},
'Si':{'Z':14,'L':3,'Zv':4,'Rv':1.97,'Rc':1.11,'m':28.09,'V':31.98,'En':2.14,'alapha':5.38,'IP':8.152,'EA':1.385},
'P':{'Z':15,'L':3,'Zv':5,'Rv':1.85,'Rc':1.06,'m':30.97,'V':26.52,'En':2.52,'alapha':3.63,'IP':10.487,'EA':0.747},
'S':{'Z':16,'L':3,'Zv':6,'Rv':1.80,'Rc':1.02,'m':32.07,'V':24.43,'En':2.96,'alapha':2.90,'IP':10.360,'EA':2.077},
'Cl':{'Z':17,'L':3,'Zv':7,'Rv':1.75,'Rc':0.99,'m':35.45,'V':22.45,'En':3.48,'alapha':2.18,'IP':12.968,'EA':3.613},
'K':{'Z':19,'L':4,'Zv':1,'Rv':2.75,'Rc':1.96,'m':39.10,'V':87.11,'En':0.45,'alapha':43.4,'IP':4.341,'EA':0.501},
'Ca':{'Z':20,'L':4,'Zv':2,'Rv':0.0,'Rc':1.74,'m':40.08,'V':0.0,'En':0.95,'alapha':22.8,'IP':6.113,'EA':0.018},
'Cr':{'Z':24,'L':4,'Zv':6,'Rv':2.20,'Rc':1.27,'m':52.00,'V':44.60,'En':1.66,'alapha':11.60,'IP':6.767,'EA':0.666},
'Mn':{'Z':25,'L':4,'Zv':7,'Rv':2.18,'Rc':1.39,'m':54.94,'V':43.40,'En':2.20,'alapha':9.40,'IP':7.434,'EA':0.0},
'Fe':{'Z':26,'L':4,'Zv':8,'Rv':2.14,'Rc':1.25,'m':55.85,'V':41.05,'En':2.20,'alapha':8.40,'IP':7.902,'EA':1.151},
'Co':{'Z':27,'L':4,'Zv':9,'Rv':2.03,'Rc':1.26,'m':58.93,'V':35.04,'En':2.56,'alapha':7.50,'IP':7.881,'EA':0.662},
'Ni':{'Z':28,'L':4,'Zv':10,'Rv':1.60,'Rc':1.21,'m':58.69,'V':17.16,'En':1.94,'alapha':6.80,'IP':7.640,'EA':1.156},
'Cu':{'Z':29,'L':4,'Zv':11,'Rv':1.40,'Rc':1.38,'m':63.55,'V':11.49,'En':1.95,'alapha':6.10,'IP':7.723,'EA':1.235},
'Zn':{'Z':30,'L':4,'Zv':12,'Rv':1.39,'Rc':1.31,'m':65.39,'V':11.25,'En':2.23,'alapha':7.10,'IP':9.394,'EA':0.0},
'Ga':{'Z':31,'L':4,'Zv':3,'Rv':1.87,'Rc':1.26,'m':69.72,'V':27.39,'En':2.42,'alapha':8.12,'IP':5.999,'EA':0.300},
'Ge':{'Z':32,'L':4,'Zv':4,'Rv':1.90,'Rc':1.22,'m':72.61,'V':28.73,'En':2.62,'alapha':6.07,'IP':7.900,'EA':1.233},
'As':{'Z':33,'L':4,'Zv':5,'Rv':1.85,'Rc':1.19,'m':74.92,'V':26.52,'En':2.82,'alapha':4.31,'IP':9.815,'EA':0.810},
'Se':{'Z':34,'L':4,'Zv':6,'Rv':1.90,'Rc':1.16,'m':78.96,'V':28.73,'En':3.01,'alapha':3.73,'IP':9.752,'EA':2.021},
'Br':{'Z':35,'L':4,'Zv':7,'Rv':1.95,'Rc':1.14,'m':79.90,'V':31.06,'En':3.22,'alapha':3.05,'IP':11.814,'EA':3.364},
'Rb':{'Z':37,'L':5,'Zv':1,'Rv':0.0,'Rc':2.11,'m':85.47,'V':0.0,'En':0.31,'alapha':47.3,'IP':4.177,'EA':0.486},
'Sr':{'Z':38,'L':5,'Zv':2,'Rv':0.0,'Rc':1.92,'m':87.62,'V':0.0,'En':0.72,'alapha':27.6,'IP':5.695,'EA':0.110},
'Mo':{'Z':42,'L':5,'Zv':6,'Rv':2.00,'Rc':1.45,'m':95.94,'V':33.51,'En':1.15,'alapha':12.80,'IP':7.092,'EA':0.746},
'Ag':{'Z':47,'L':5,'Zv':11,'Rv':1.72,'Rc':1.53,'m':107.87,'V':21.31,'En':1.83,'alapha':7.20,'IP':7.576,'EA':1.302},
'Cd':{'Z':48,'L':5,'Zv':12,'Rv':1.58,'Rc':1.48,'m':112.41,'V':16.52,'En':1.98,'alapha':7.20,'IP':8.994,'EA':0.0},
'In':{'Z':49,'L':5,'Zv':3,'Rv':1.93,'Rc':1.44,'m':114.82,'V':30.11,'En':2.14,'alapha':10.20,'IP':5.786,'EA':0.300},
'Sn':{'Z':50,'L':5,'Zv':4,'Rv':2.22,'Rc':1.41,'m':118.71,'V':45.83,'En':2.30,'alapha':7.70,'IP':7.344,'EA':1.112},
'Sb':{'Z':51,'L':5,'Zv':5,'Rv':2.10,'Rc':1.38,'m':121.76,'V':38.79,'En':2.46,'alapha':6.60,'IP':8.64,'EA':1.07},
'Te':{'Z':52,'L':5,'Zv':6,'Rv':2.06,'Rc':1.35,'m':127.60,'V':36.62,'En':2.62,'alapha':5.50,'IP':9.01,'EA':1.971},
'I':{'Z':53,'L':5,'Zv':7,'Rv':2.10,'Rc':1.33,'m':126.90,'V':38.79,'En':2.78,'alapha':5.35,'IP':10.451,'EA':3.059},
'Gd':{'Z':64,'L':6,'Zv':10,'Rv':2.59,'Rc':1.79,'m':157.25,'V':72.78,'En':2.00,'alapha':23.50,'IP':6.15,'EA':0.50},
'Pt':{'Z':78,'L':6,'Zv':10,'Rv':1.75,'Rc':1.28,'m':195.08,'V':22.45,'En':2.28,'alapha':6.50,'IP':9.00,'EA':2.128},
'Au':{'Z':79,'L':6,'Zv':11,'Rv':1.66,'Rc':1.44,'m':196.97,'V':19.16,'En':2.65,'alapha':5.80,'IP':9.226,'EA':2.309},
'Hg':{'Z':80,'L':6,'Zv':12,'Rv':1.55,'Rc':1.49,'m':200.59,'V':15.60,'En':2.20,'alapha':5.70,'IP':10.438,'EA':0.0},
'Tl':{'Z':81,'L':6,'Zv':3,'Rv':1.96,'Rc':1.48,'m':204.38,'V':31.54,'En':2.25,'alapha':7.60,'IP':6.108,'EA':0.200},
'Pb':{'Z':82,'L':6,'Zv':4,'Rv':2.02,'Rc':1.47,'m':207.20,'V':34.53,'En':2.29,'alapha':6.80,'IP':7.417,'EA':0.364},
'Bi':{'Z':83,'L':6,'Zv':5,'Rv':2.10,'Rc':1.46,'m':208.98,'V':38.79,'En':2.34,'alapha':7.40,'IP':7.289,'EA':0.946},

}



def GetAbsoluteAtomicProperty(element='C',propertyname='m'):
    """
    Get the absolute property value with propertyname for the given atom.
    """
    
    PropertyDic=AtomProperty[element]
    
    return PropertyDic[propertyname]


def GetRelativeAtomicProperty(element='C',propertyname='m'):
    """
    Get the absolute property value with propertyname for the given atom.
    """
    
    CpropertyDic=AtomProperty['C']
    PropertyDic=AtomProperty[element]
    
    return PropertyDic[propertyname]/(CpropertyDic[propertyname]+0.0)

###############################################################################

if __name__=="__main__":
    
    for i,j in AtomProperty.items():
        print j
    print GetAbsoluteAtomicProperty(element='S',propertyname='En')
    print GetRelativeAtomicProperty(element='S',propertyname='En')
