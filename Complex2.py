from cmath import sqrt
import math
import Complex as Bas;

def SumV(Data1,Data2):
    res=[]
    if len(Data1)==len(Data2):
        for i in range(len(Data1)):
            res.append(Bas.sum(Data1[i][0],Data1[i][1],Data2[i][0],Data2[i][1]))
        return res
    else: return res
def InvV(Data1):
    res=[]
    for i in range(len(Data1)):
        res.append(Bas.pro(Data1[i][0],Data1[i][1],-1,0))
    return res
def eVec(Escalar,Data1):
    res=[]
    for i in range(len(Data1)):
        res.append(Bas.pro(Data1[i][0],Data1[i][1],Escalar[0],Escalar[1]))
    return res
def SumM(Data1,Data2):
    res=[]
    if len(Data1)==len(Data2) and len(Data1[0])==len(Data2[0]):
        for i in range(len(Data1)):
            columna=[]
            for j in range(len(Data1[i])):
                columna.append(Bas.sum(Data1[i][j][0],Data1[i][j][1],Data2[i][j][0],Data2[i][j][1]))
            res.append(columna)
        return res
    else: return res
def InvM(Data1):
    res=[]
    for i in range(len(Data1)):
        columna=[]
        for j in range(len(Data1[i])):
            columna.append(Bas.pro(Data1[i][j][0],Data1[i][j][1],-1,0))
        res.append(columna)
    return res
def eMat(Escalar,Data1):
    res=[]
    for i in range(len(Data1)):
        columna=[]
        for j in range(len(Data1[i])):
            columna.append(Bas.pro(Data1[i][j][0],Data1[i][j][1],Escalar[0],Escalar[1]))
        res.append(columna)
    return res
def Tra(Data1):
    res=[]
    for i in range(len(Data1)):
        columna=[]
        for j in range(len(Data1[i])):
            columna.append(Bas.pro(Data1[j][i][0],Data1[j][i][1],1,0))
        res.append(columna)
    return res
def Conj(Data1):
    res=[]
    for i in range(len(Data1)):
        columna=[]
        for j in range(len(Data1[i])):
            columna.append(Bas.conj(Data1[i][j][0],Data1[i][j][1]))
        res.append(columna)
    return res
def Adj(Data1):
    res=[]
    for i in range(len(Data1)):
        columna=[]
        for j in range(len(Data1[i])):
            columna.append(Bas.conj(Data1[j][i][0],Data1[j][i][1]))
        res.append(columna)
    return res
def MulM(Data1,Data2):
    res=[]
    for i in range(len(Data1)):
        columna=[]
        for j in range(len(Data1[i])):
            Rp=0
            Ip=0
            for k in range(len(Data2)):
                Rp+=Data1[i][k][0]*Data2[k][j][0]-Data1[i][k][1]*Data2[k][j][1]
                Ip+=Data1[i][k][0]*Data2[k][j][1]+Data2[k][j][0]*Data1[i][k][1]
            Num=(round(Rp,2),round(Ip,2))
            columna.append(Bas.pro(Num[0],Num[1],1,0))
        res.append(columna)
    return res
def AcMV(Data1,Data2):
    res=[]
    for i in range(len(Data1)):
        columna=[]
        Rp=0
        Ip=0
        for k in range(len(Data2)):
            Rp+=Data1[i][k][0]*Data2[k][0]-Data1[i][k][1]*Data2[k][1]
            Ip+=Data1[i][k][0]*Data2[k][1]+Data2[k][0]*Data1[i][k][1]
        Num=(Rp,Ip)
        columna.append(Bas.pro(Num[0],Num[1],1,0))
        res.append(columna)
    return res
def PriV(Data1,Data2):
    res=0
    Rp=0
    Ip=0
    for i in range(len(Data1)):
        Rp+=Data1[i][0]*Data2[i][0]-((-1)*Data1[i][1])*Data2[i][1]
        Ip+=Data1[i][0]*Data2[i][1]+Data2[i][0]*((-1)*Data1[i][1])
    res=Bas.pro(Rp,Ip,1,0)
    return res
def NorV(Data1):
    res=0
    Rp=0
    for i in range(len(Data1)):
        Rp+=Data1[i][0]**2+Data1[i][1]**2
    res=Rp
    return round(math.sqrt(res),2)
def DisV(Data1,Data2):
    res=[]
    for i in range(len(Data1)):
        res.append((Data1[i][0]-Data2[i][0],Data1[i][1]-Data2[i][1]))
    d=0
    for j in range(len(res)):
        d+=res[j][0]**2+res[j][1]**2
    return round(math.sqrt(d),2)
def Unit(Data1):
    res=[]
    Un=[]
    for i in range(len(Data1)):
        colum=[]
        col1=[]
        for j in range (len(Data1[i])):
            if(i==j):col1.append(1)
            else: col1.append(0)
            colum.append((Data1[j][i][0],Data1[j][i][1]*-1))
        Un.append(col1)
        res.append(colum)
    res=MulM(res,Data1)
    if(res==Un):return True
    else:return False
def Herm(Data1):
    res=[]
    Un=[]
    for i in range(len(Data1)):
        colum=[]
        col1=[]
        for j in range (len(Data1[i])):
            if(i==j):col1.append(1)
            else: col1.append(0)
            colum.append((Data1[j][i][0],Data1[j][i][1]*-1))
        Un.append(col1)
        res.append(colum)
    if(res==Data1):return True
    else:return False
def Prot(Data1,Data2):
    res=[]
    if(type(Data1[0])==list):
        for i in Data1:
            for j in i:
                if (type(Data2[0])==list):
                    res.append(eMat(j,Data2))
                else:
                    res.append(eVec(j,Data2))
    else:
        for i in Data1:
            if (type(Data2[0])==list):
                res.append(eMat(i,Data2))
            else:
                res.append(eVec(i,Data2))
    return res

print(SumV([(1,2),(2,2),(3,3)],[(3,2),(1,1),(1,1)]))
print(InvV([(1,2),(2,2),(3,3),(3,2),(1,1),(1,1)]))
print(eVec((1,2),[(1,2),(2,2),(3,3),(3,2),(1,1),(1,1)]))
print(SumM([[(1,2),(2,1)],[(3,2),(5,4)]],[[(1,2),(2,1)],[(3,2),(5,4)]]))
print(InvM([[(1,2),(2,2)],[(3,3),(3,2)],[(1,1),(1,1)]]))
print(eMat((1,2),[[(1,2),(2,2)],[(3,3),(3,2)],[(1,1),(1,1)]]))
print(Tra([[(1,2),(2,2),(5,5)],[(3,3),(3,2),(4,4)],[(1,1),(1,1),(6,4)]]))
print(Conj([[(1,2),(2,2),(5,5)],[(3,3),(3,2),(4,4)],[(1,1),(1,1),(6,4)]]))
print(Adj([[(1,2),(2,2),(5,5)],[(3,3),(3,2),(4,4)],[(1,1),(1,1),(6,4)]]))
print(MulM([[(1,1),(2,2),(3,3)],[(3,3),(4,4),(5,5)],[(3,3),(4,4),(5,5)]],[[(1,1),(2,2),(6,6)],[(3,3),(4,4),(5,5)],[(3,3),(4,4),(5,5)]]))
print(AcMV([[(1,1),(2,2),(3,3)],[(3,3),(4,4),(5,5)],[(3,3),(4,4),(5,5)]],[(1,1),(2,2),(3,3)]))
print(PriV([(1,2),(2,2),(3,3)],[(3,2),(1,1),(1,1)]))
print(NorV([(4,3),(6,-4),(14,-7),(0,13)]))
print(DisV([(3,0),(1,0),(2,0)],[(2,0),(2,0),(-1,0)]))
print(Unit([[(1,0),(0,0)],[(0,0),(1,0)]]))
print(Herm([[(1,0),(0,-1)],[(0,1),(1,0)]]))
print(Prot([[(1,1),(2,2),(3,3)]],[[(1,1),(2,2),(6,6)],[(3,3),(4,4),(5,5)],[(3,3),(4,4),(5,5)]]))