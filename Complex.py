import math

def sum(Real_P1,Img_P1,Real_P2,Img_P2):
    res= (Real_P1+Real_P2,Img_P1+Img_P2)
    if (res[1]==0): return res[0]
    elif res[1]>0: return "{}+{}i".format(res[0],res[1])
    else: return "{}{}i".format(res[0],res[1])
def res(Real_P1,Img_P1,Real_P2,Img_P2):
    res= (Real_P1-Real_P2,Img_P1-Img_P2)
    if (res[1]==0): return res[0]
    elif res[1]>0: return "{}+{}i".format(res[0],res[1])
    else: return "{}{}i".format(res[0],res[1])
def pro(Real_P1,Img_P1,Real_P2,Img_P2):
    res= (Real_P1*Real_P2-Img_P1*Img_P2,Real_P1*Img_P2+Real_P2*Img_P1)
    if (res[1]==0): return res[0]
    elif res[1]>0: return "{}+{}i".format(res[0],res[1])
    else: return "{}{}i".format(res[0],res[1])
def div(Real_P1,Img_P1,Real_P2,Img_P2):
    res= (round((((Real_P1*Real_P2)+(Img_P1*Img_P2))/(Real_P2**2+Img_P2**2)),2),round((((Real_P2*Img_P1)-(Real_P1*Img_P2))/(Real_P2**2+Img_P2**2)),2))
    if (res[1]==0): return res[0]
    elif res[1]>0: return "{}+{}i".format(res[0],res[1])
    else: return "{}{}i".format(res[0],res[1])
def mod(Real_P1,Img_P1):
    res=(round(math.sqrt(Real_P1**2+Img_P1**2),2))
    return res
def conj(Real_P1,Img_P1):
    res=(Real_P1,(-1)*Img_P1)
    if (res[1]==0): return res[0]
    elif res[1]>0: return "{}+{}i".format(res[0],res[1])
    else: return "{}{}i".format(res[0],res[1])
def ctop(Real_P1,Img_P1):
    res=(round(math.sqrt(Real_P1**2+Img_P1**2),2),round(math.atan(Img_P1/Real_P1)*180/math.pi,2))
    return "({},{}°)".format(res[0],res[1])
def fase(Real_P1,Img_P1):
    res= round(math.atan(Img_P1/Real_P1)*180/math.pi,2)
    return "{}°".format(res)
