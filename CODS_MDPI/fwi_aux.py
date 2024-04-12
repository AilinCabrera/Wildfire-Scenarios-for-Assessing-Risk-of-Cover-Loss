import os
import locale
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
locale.setlocale(locale.LC_TIME, '')
tipo_data = 1
#################################################
def FFMC(T, H, W, ro, FO):
    mo = (147.2 * (101 - FO))/(59.5 + FO)
    
    if ro > 0.5:
        rf = ro - 0.5
        
        if mo <= 250:
            mr = mo + 42.5 * rf * np.exp(-100/(251 - mo)) * (1 - np.exp(-6.93/rf))
        
        if mo > 250:
            mr = mo + 42.5 * rf * np.exp(-100/(251 - mo)) * (1 - np.exp(-6.93/rf)) + 0.0015 * (mo - 150)**2 * rf**0.5
            
        mo = mr
    
    Ed = 0.942 * H**0.679 + 11 * np.exp((H - 100)/10) + 0.18 * (21.1 - T) * (1 - np.exp(-0.115 * H))

    if mo > Ed:
        ko = 0.424 * (1 - (H/100)**1.7) + 0.0694 * W**0.5 *  (1 - (H/100)**8)
        kd = ko * 0.581 * np.exp(0.0365 * T)
        m  = Ed + (mo - Ed) * 10**-kd

    if mo < Ed:
        Ew = 0.618 * H**0.753 + 10 * np.exp((H - 100)/10) + 0.18 * (21.1 - T) * (1 - np.exp(-0.115 * H))
        
        if mo < Ew:
            kl = 0.424 * (1 - ((100 - H)/100)**1.7) + 0.0694 * W**0.5 * (1 - ((100 - H)/100)**8)
            kw = kl * 0.581 * np.exp(0.0365 * T)
            m  = Ew - (Ew - mo) * 10**-kw
            
        if Ed >= mo and mo >= Ew:
            m = mo

    F  = 59.5 * (250 - m)/(147.2 + m)
    return(m,np.round(F,1))

#################################################
def DMC(T, H, Le, ro, po):
    if ro > 1.5:
        re = 0.92 * ro - 1.27
        
        mo = 20 + np.exp(5.6348 - (po/43.43))
        
        if po <= 33:
            b = 100/(0.5+ 0.3*po)
        
        if po > 33 and po <= 65:
            b = 14 - 1.3 * np.log(po)
            
        if po > 65:
            b = 6.2 * np.log(po) - 17.2                
            
        Mr = mo + (1000 * re/(48.77 + b * re))         
        Pr = 244.72 - 43.43 * np.log(Mr - 20)          
        po = Pr                                        
        
    K = 1.894 * (T + 1.1) * (100 - H) * (Le  *10**-6)
    P = po + 100 * K                             
    return(P)

#################################################
def filtrar_mes_DCM(mes,i):
    return(9)

#################################################
def DC(T, H, ro, Do, Lf):
  if ro > 2.8:                              # CAMBIO           r
     rd = 0.83*ro - 1.27
     Qo = 800 * np.exp(-Do/400)             #CAMBIO             Qo = 800 * np.exp(Do/400)
     Qr = Qo + 3.937 * rd
     Dr = 400 * np.log(800/Qr)
     Do = Dr                                #Agregar
       
  V = 0.36 * (T + 2.8) + Lf
  D = Do + 0.5 * V
  return(np.abs(D))

#################################################
def filtrar_mes_lf(mes,i):
    return(1.4)

def ISI(W,m):
  fw = np.exp(0.05039 * W)
  fF = (91.9 * np.exp(-0.1386 * m)) * (1 + m**5.31/(4.93*10**7))   #cambio fF = (91.9 * np.exp(-0.1386)) * ((1+ m *exp(5.31))/(4.93*107))
  R = 0.208 * fw * fF
  return(R)


def BUI(DMC, DC):
    if DMC <= 0.4 * DC:                                       
        U = 0.8 * DMC * DC / (DMC + 0.4* DC)                                   
    else:
        U = DMC - (1 - 0.8 * DC/(DMC + 0.4  * DC ))* (0.92 + (0.0114 * DMC)**1.7) 
    if U <=0:
        U = 0
    return(np.round(U,1))


def FWI2(BUI,SI):
    if BUI >= 80:
        BB = 0.1 * SI * (1000/(25+108.64/np.exp(0.023 * BUI)))
    else:
        BB = 0.1 * SI * (0.626 * BUI**0.809 + 2)

    if BB-1 <= 0:
        FWI = BB
    else:
        SL = 2.72*(0.434 * np.log(BB))**0.647
        FWI = np.exp(SL)

    return(FWI)

def DSR(S):
  DSR = 0.0272 * (S**1.77)           #CAMBIO             DSR = 0.0272 * (FWI() ^ 1.77) 
  return(DSR)



def xxx(data_csv):
    T   = data_csv['TEMP'] # °C 
    H   = data_csv['RH']   # %
    W   = data_csv['WIND'] # km/h
    ro  = data_csv['RAIN'] # mm
    mes = data_csv['MES']
    
    if tipo_data == 1:
          cont     = 0
          FFMC_LIS = []
          M_LIS    = []
          for i in np.arange(0,len(data_csv),1):
              if cont == 0:
                  FO = 85

              m,FFMC_OUT = FFMC(T = T[i], H = H[i], W = W[i], ro = ro[i], FO = FO)
              FO         = FFMC_OUT ; cont += 1
              M_LIS      = np.append(M_LIS,m)
              FFMC_LIS   = np.append(FFMC_LIS,FFMC_OUT)

          dff = pd.DataFrame(FFMC_LIS, columns = ['FFMC'])
          dfm = pd.DataFrame(M_LIS, columns = ['M'])
          df  = pd.merge(data_csv.reset_index(),dfm.reset_index(), left_index=False, right_index=False) ; df.pop('index')
          df  = pd.merge(df.reset_index(),dff.reset_index(), left_index=False, right_index=False) ; df.pop('index')


    if tipo_data == 1:
        cont     = 0
        DCM_LIS = []
        for i in np.arange(0,len(data_csv),1):
            if cont == 0:
                po = 6
                
            Le = filtrar_mes_DCM(mes,i)
            DCM_OUT = DMC(T = T[i], H = H[i],Le = Le ,ro = ro[i], po = po)
            po = DCM_OUT ; cont += 1
            DCM_LIS = np.append(DCM_LIS,np.round(DCM_OUT,1))
            
        df2 = pd.DataFrame(DCM_LIS, columns = ['DCM'])
        df = pd.merge(df.reset_index(),df2.reset_index(), left_index=False, right_index=False) ; df.pop('index')


    if tipo_data == 1:
        cont     = 0
        DC_LIS = []
        for i in np.arange(0,len(data_csv),1):
            if cont == 0:
                Do = 15
                
            Lf = filtrar_mes_lf(mes,i)
            DC_OUT = DC(T = T[i], H = H[i], ro = ro[i], Do = Do, Lf = Lf)
            Do = DC_OUT ; cont += 1
            DC_LIS = np.append(DC_LIS,np.round(DC_OUT,1))
        
        df2 = pd.DataFrame(DC_LIS, columns = ['DC'])
        df = pd.merge(df.reset_index(),df2.reset_index(), left_index=False, right_index=False) ; df.pop('index')


    if tipo_data == 1:
        ISI_LIS = []
        for i in np.arange(0,len(data_csv),1):

            ISI_OUT = ISI(W = W[i], m = df['M'][i])
            ISI_LIS = np.append(ISI_LIS,np.round(ISI_OUT,1))
        
        df2 = pd.DataFrame(ISI_LIS, columns = ['ISI'])
        df = pd.merge(df.reset_index(),df2.reset_index(), left_index=False, right_index=False) ; df.pop('index')


    if tipo_data == 1:
        BUI_LIS = []
        for i in np.arange(0,len(data_csv),1):

            BUI_OUT = BUI(DMC = df['DCM'][i],DC = df['DC'][i])
            BUI_LIS = np.append(BUI_LIS,np.round(BUI_OUT,1))
        
        df2 = pd.DataFrame(BUI_LIS, columns = ['BUI'])
        df = pd.merge(df.reset_index(),df2.reset_index(), left_index=False, right_index=False) ; df.pop('index')


    df3 = df
    if tipo_data == 1:
        FWI_LIS = []
        for i in np.arange(0,len(data_csv),1):

            FWI_OUT = FWI2(BUI = df['BUI'][i], SI = df['ISI'][i])
            FWI_LIS = np.append(FWI_LIS,np.round(FWI_OUT,1))
        
        df2 = pd.DataFrame(FWI_LIS, columns = ['FWI'])
        df3 = pd.merge(df3.reset_index(),df2.reset_index(), left_index=False, right_index=False) ; df3.pop('index')


    df4 = df3
    if tipo_data == 1:
        DSR_LIS = []
        for i in np.arange(0,len(data_csv),1):
            DSR_OUT = DSR(S = df3['FWI'][i])
            DSR_LIS = np.append(DSR_LIS,np.round(DSR_OUT,2))
        
        df2 = pd.DataFrame(DSR_LIS, columns = ['DSR'])
        df4 = pd.merge(df4.reset_index(),df2.reset_index(), left_index=False, right_index=False) ; df4.pop('index')


    df4 = df4.interpolate()


    return df4