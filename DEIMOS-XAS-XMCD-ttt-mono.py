# -*- coding: utf-8 -*-
"""
FR
Permet de normaliser, recentrer sur 0, combiner et les spectres XAS/XMCD issues
de DEIMOS-SOLEIL puis d'appliquer automatiquement les règles de somme.
Un numéro de scan en entrée (avec un paramètre Tz sans incidence sur les calculs),
les 7 suivants étant le même échantillon avec +/-H et circulaire droite (CR) et gauche (CL).

EN
Allow the import, normalisation, centring and merging of XAS/XMCD spectra from
DEIMOS-SOLEIL beamline to automatically apply sum rules.
One scan number to give (with a Tz parameter without incidence on calculations),
the next 7 should be the same sample with +/-H and circular right (CR) and left (CL).



TO DO: use numpy.array instead of lists
"""
import numpy as np
import sys
#import matplotlib
import matplotlib.pyplot as plt


""" --------------- running mean function ---------------------- """
def running_mean(x, i, N):
    out_mean=0
    for j in range(-N,N+1):
        out_mean=out_mean+x[j+i]
    out_mean=(out_mean-x[i])/(2*N)
    return out_mean
""" ----------------------------------------------------------- """


""" -------- function find nearest value in array --------------- """
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
""" ----------------------------------------------------------- """



""" --------------------- parameters ---------------------------- """
DebutScan=75
Tz=27

MethodePH=False
NbScans=8
""" ----------------------------------------------------------- """


""" ---------------- advanced parameters ----------------------- """
ScansPos=[0,3,5,6]
ScansNeg=[1,2,4,7]
NbScans=len(ScansPos)+len(ScansNeg)
""" ----------------------------------------------------------- """





"""--------------- error checking ---------------"""
if len(ScansPos) != len(ScansNeg): # useless for now
    sys.exit('Error: Number of positive and negative scans is not equal')
""" ----------------------------------------------------------- """




"""--------------- import of data ---------------"""
energy=[]
i0=[]
it=[]
itio=[]     # TEY
field=[]
if1=[]
ifio=[]     # fluo
ition=[]
ifion=[]

minEnergy=10000
maxEnergy=0

for file in range(0,8):
    energy.append([])
    i0.append([])
    it.append([])
    itio.append([])
    field.append([])
    if1.append([])
    ifio.append([])
    ition.append([])
    ifion.append([])

    fichier = 'scan_'+ '%03i'%(DebutScan+file)+'.txt'
    data = np.genfromtxt('./2020-02-12/'+fichier, skip_header=1) #delimiter='        '
    for i in range(0,len(data)):
        energy[file].append(data[i,2])
        i0[file].append(data[i,6])
        it[file].append(data[i,7])
        itio[file].append(data[i,8])
        field[file].append(data[i,10])
        if1[file].append(data[i,12])
        ifio[file].append(data[i,13])
    if min(energy[file])<minEnergy:
        minEnergy=min(energy[file])
    if max(energy[file])>maxEnergy:
        maxEnergy=max(energy[file])
NbPts=len(data)
""" ----------------------------------------------------------- """






"""--------------------- edge detection ---------------------"""
Seuil='unknown edge'
nh3d=0
if abs(energy[0][0]-705)<20: # Fe
    Seuil='Fe'
    nh3d=3.7  # 3.7 for L10, 3.39 for bulk
    picL3tab=705
    picL2tab=718

elif abs(energy[0][0]-775)<20: # Co
    Seuil='Co'
    nh3d=2.49
    picL3tab=775
    picL2tab=790
""" ----------------------------------------------------------- """





"""--------------------- merge + normalisation ---------------------"""
energyCR=[]
energyCL=[]
energyCLR=[]
XASTEY=[]
XMCDTEY=[]
XASFluo=[]
XMCDFluo=[]
XASTEYn=[]  # normalised
XMCDTEYn=[]  # normalised
XASFluon=[]  # normalised
XMCDFluon=[]  # normalised

mu_R_TEY=[]
mu_L_TEY=[]
mu_R_Fluo=[]
mu_L_Fluo=[]
mu_R_TEYn=[] # normalised
mu_L_TEYn=[] # normalised
mu_R_Fluon=[] # normalised
mu_L_Fluon=[] # normalised

for i in range(0,NbPts):
    energyCR.append((energy[0][i]+energy[3][i]+energy[5][i]+energy[6][i])/4) # 0123           4567
    energyCL.append((energy[1][i]+energy[2][i]+energy[4][i]+energy[7][i])/4) # +--+ and -H so -++- for Fe
    energyCLR.append((energyCR[i]+energyCL[i])/2)
    
for file in range(0,8):
    normDebutTEY=running_mean(itio[file],10,5)
    normFinTEY=running_mean(itio[file],-10,5)
    normDebutFluo=running_mean(ifio[file],10,5)
    normFinFluo=running_mean(ifio[file],-10,5)
    DeltaNormTEY=normFinTEY-normDebutTEY
    DeltaNormFluo=normFinFluo-normDebutFluo
    if (normDebutTEY>normFinTEY or normDebutFluo>normFinFluo) and MethodePH==True:
        DeltaNormTEY=1
        DeltaNormFluo=1
    if (normDebutTEY>normFinTEY or normDebutFluo>normFinFluo) and MethodePH==False:
        DeltaNormTEY=max(itio[file])
        DeltaNormFluo=max(ifio[file])
    for i in range(0,NbPts):
        ition[file].append((itio[file][i]-normDebutTEY)/DeltaNormTEY)
        ifion[file].append((ifio[file][i]-normDebutFluo)/DeltaNormFluo)

for i in range(0,NbPts): # merge
    mu_R_TEYn.append((ition[0][i]+ition[3][i]+ition[5][i]+ition[6][i])/4)
    mu_L_TEYn.append((ition[1][i]+ition[2][i]+ition[4][i]+ition[7][i])/4)
    mu_R_Fluon.append((ifion[0][i]+ifion[3][i]+ifion[5][i]+ifion[6][i])/4)
    mu_L_Fluon.append((ifion[1][i]+ifion[2][i]+ifion[4][i]+ifion[7][i])/4)
    mu_R_TEY.append((itio[0][i]+itio[3][i]+itio[5][i]+itio[6][i])/4)
    mu_L_TEY.append((itio[1][i]+itio[2][i]+itio[4][i]+itio[7][i])/4)
    mu_R_Fluo.append((ifio[0][i]+ifio[3][i]+ifio[5][i]+ifio[6][i])/4)
    mu_L_Fluo.append((ifio[1][i]+ifio[2][i]+ifio[4][i]+ifio[7][i])/4)

for i in range(0,NbPts):
    XASTEY.append((mu_R_TEY[i]+mu_L_TEY[i])/2)
    XMCDTEY.append(mu_R_TEY[i]-mu_L_TEY[i])
    XASFluo.append((mu_R_Fluo[i]+mu_L_Fluo[i])/2)
    XMCDFluo.append(mu_R_Fluo[i]-mu_L_Fluo[i])
    
    XASTEYn.append((mu_R_TEYn[i]+mu_L_TEYn[i])/2)
    XMCDTEYn.append(mu_R_TEYn[i]-mu_L_TEYn[i])
    XASFluon.append((mu_R_Fluon[i]+mu_L_Fluon[i])/2)
    XMCDFluon.append(mu_R_Fluon[i]-mu_L_Fluon[i])

if Seuil=='Fe':  # because -++- for Fe (see first loop of this part)
    for i in range(0,NbPts):
        XMCDTEYn[i]=-XMCDTEYn[i]
        XMCDFluon[i]=-XMCDFluon[i]
        XMCDTEY[i]=-XMCDTEY[i]
        XMCDFluo[i]=-XMCDFluo[i]
""" ----------------------------------------------------------- """







"""------- search max, min, calculation of step and holes for TEY ------"""
x=np.array(XASTEYn)

picL3tabnearest=find_nearest(energyCLR, picL3tab)
picL3tabi=np.where(energyCLR == picL3tabnearest)
picL3tabi=int(picL3tabi[0])

picL2tabnearest=find_nearest(energyCLR, picL2tab)
picL2tabi=np.where(energyCLR == picL2tabnearest)
picL2tabi=int(picL2tabi[0])


picL3=int(np.where(x == max(x[picL3tabi-50:picL3tabi+50]))[0])
picL2=int(np.where(x == max(x[picL2tabi-50:picL2tabi+50]))[0])


annotx=max(energyCLR) # for text position on plot
annoty=max(x)

minXASTEY=picL3
for i in range(picL3,picL2):
    if XASTEYn[i]<XASTEYn[minXASTEY]:
        minXASTEY=i

step=[]
holes=[]
for i in range(0,NbPts):
    step.append((((2/(3*np.pi)*np.arctan((energyCLR[i]-energyCLR[picL3])/0.5))+1/3)+((1/(3*np.pi)*np.arctan((energyCLR[i]-energyCLR[picL2])/0.5))+1/6))*1.007-0.01)
    holes.append(x[i]-step[i])
""" ----------------------------------------------------------- """






"""--------------------- plots ---------------------"""
plt.title('XAS, XMCD, step TEY')
plt.plot(energyCLR, x)
plt.plot(energyCLR, XMCDTEYn)
plt.plot(energyCLR, step)
plt.plot(energyCLR[picL3], x[picL3], 'x',color='black')
plt.plot(energyCLR[picL2], x[picL2], 'x',color='black')
plt.plot(energyCLR[minXASTEY], x[minXASTEY], 'x',color='red')
plt.annotate('TEY', xy=(annotx-15,annoty-2), color="black")
plt.annotate('Tz='+str(Tz), xy=(annotx-15,annoty-2.4), color="black")
plt.annotate(Seuil, xy=(annotx-15,annoty-2.8), color="black")
plt.grid(True)
plt.show()


plt.title('holes TEY')
plt.plot(energyCLR, holes)
#plt.plot(energyCLR, XASTEYn)
plt.plot(energyCLR[picL3], holes[picL3], 'x',color='black')
plt.plot(energyCLR[picL2], holes[picL2], 'x',color='black')
plt.plot(energyCLR[minXASTEY], holes[minXASTEY], 'x',color='red')
plt.annotate('TEY', xy=(annotx-15,annoty-2), color="black")
plt.annotate('Tz='+str(Tz), xy=(annotx-15,annoty-2.4), color="black")
plt.annotate(Seuil, xy=(annotx-15,annoty-2.8), color="black")
plt.grid(True)
plt.show()
""" ----------------------------------------------------------- """






"""--------------------- integration and sum rules ---------------------"""
IL3 = np.trapz(holes[0:minXASTEY],energyCLR[0:minXASTEY], dx=1)
IL2 = np.trapz(holes[minXASTEY:-1],energyCLR[minXASTEY:-1], dx=1)
DeltaIL3 = np.trapz(XMCDTEYn[0:minXASTEY],energyCLR[0:minXASTEY], dx=1)
DeltaIL2 = np.trapz(XMCDTEYn[minXASTEY:-1],energyCLR[minXASTEY:-1], dx=1)
C=nh3d/(IL3+IL2)
m_orb_d=-2*C*(DeltaIL3+DeltaIL2)/3
m_spineff_d=-C*(DeltaIL3-2*DeltaIL2)

print(' Tz='+str(Tz)+'\n Edge='+Seuil+'\n L3peak= '+str(energyCLR[picL3])+'\n L2peak='+str(energyCLR[picL2])+'\n min='+str(energyCLR[minXASTEY])+'\n nh3d='+str(nh3d)+'\n C='+str(C)+'\n morb='+str(m_orb_d)+'\n mspin eff='+str(m_spineff_d)+'\n mspin eff/morb='+str(m_spineff_d/m_orb_d))
print('-----------------------------------------')

print('----- to copy and paste directly in Origin workbook -----')
print(str(Tz)+'\t'+Seuil+'\t'+str(energyCLR[picL3])+'\t'+str(energyCLR[picL2])+'\t'+str(energyCLR[minXASTEY])+'\t'+str(nh3d)+'\t'+str(C)+'\t'+str(m_orb_d)+'\t'+str(m_spineff_d)+'\t'+str(m_spineff_d/m_orb_d))
print('-----------------------------------------')

""" ----------------------------------------------------------- """

# Copy and paste for the column title in Origin Workbook:
# Tz	Edge	L3 peak	L2 peak	interpeak min	nh3d	C	morb	mspin eff	mspin eff/morb



"""
write .txt file as output
"""
#f=open('./output/'+Seuil+'-' +'Tz'+'%3.1f'%Tz+ '.txt','w+')
#f.write('E \t XAS TEY \t XMCD TEY \t XAS fluo \t XMCD fluo \t step TEY \t holes TEY \n')
#f.write('  \t   \t   \t   \t   \t   \t   \n')
#f.write('E ('+Seuil+' edge) \t Tz='+str(Tz)+' \t Tz='+str(Tz)+' \t Tz='+str(Tz)+' \t Tz='+str(Tz)+' \t Tz='+str(Tz)+' \t Tz='+str(Tz)+' \n')
#for i in range(0,NbPts):
#    f.write('%.6E \t %.6E \t %.6E \t %.6E \t %.6E \t %.6E \t %.6E \n' %(energyCLR[i], XASTEYn[i], XMCDTEYn[i], XASFluon[i], XMCDFluon[i], step[i], holes[i]))
#f.close()
#
#
# f=open('./output/non-norme-'+Seuil+'-' +'Tz'+'%3.1f'%Tz+ '.txt','w+')
# f.write('E \t XAS TEY \t XMCD TEY \t XAS fluo \t XMCD fluo \n')
# f.write('  \t   \t   \t   \t   \n')
# f.write('E ('+Seuil+' edge) \t Tz='+str(Tz)+' \t Tz='+str(Tz)+' \t Tz='+str(Tz)+' \t Tz='+str(Tz)+' \n')
# for i in range(0,NbPts):
#     f.write('%.6E \t %.6E \t %.6E \t %.6E \t %.6E \t %.6E \n' %(energyCLR[i], step[i], XASTEY[i], -XMCDTEY[i], XASFluo[i], -XMCDFluo[i]))
# f.close()




print('done')





