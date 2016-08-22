import random as ran
from scipy import interpolate
#import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline 
ran.seed(1) # Set the random nymber generator to a fixed sequence.


def ar1(langd,phi,c=0):
    lista=[]
    x=0.0
    for i in range(langd):
        x=c+phi*x + ran.gauss(0,1)
        lista.append(x)
    return lista

def binom(langd,chance,degeneration=1.0):
    lista=[0.0]
    for i in range(1,langd+1):
        x=lista[-1] + (ran.uniform(0,1)<chance*(degeneration**i))
        lista.append(x)
    return lista[1:]

def binom2(langd,chance,degeneration=1.0):
    lista=binom(langd,chance,degeneration)
    listb=[]
    for i in range(1,langd+1):
        x=lista[i-1]
        listb.append(x/float(i))
    return listb

def choose(n,k):
    k=int(k)
    n=int(n)
    t=1
    if k>n-k:
        for j in range(k+1,n+1):
            t*=j
        for j in range(1,n-k+1):
            t/=j
    else:
        for j in range(n-k+1,n+1):
            t*=j
        for j in range(1,k+1):
            t/=j
    return t


class rateDistribution:
    def __init__(self,rank,antal):
        self.rank=float(rank)
        self.antal=float(antal)
    def cumulativeDensity(self,y,expectedvalue=0):
        konstant=self.rank*choose(self.antal,self.rank)
        summa=0
        for i in range(int(self.antal-self.rank+1)):
            summa+=((-1)**i)*choose(self.antal-self.rank,i)*(y**(i+self.rank+expectedvalue))/float(i+self.rank+expectedvalue)
        return konstant*summa
    def expectedValue(self):
        return self.rank/(self.antal+1)
    def impOrder(self):
        return self.rank+1



kons=1.0
def multiplication(lista):
    """bygger pa ett enkelt matematisk samband"""
    rank=0
    antal=0
    for i in lista:
        rank*=kons
        antal*=kons
        rank+=i.rank
        antal+=i.antal
    rd=rateDistribution(rank,antal)
    return rd
            

class dataStruct:
    def __init__(self,nr,langd):
        self.nr=nr
        self.langd=langd
        self.matris=[]
        self.start("binom2",[self.langd,0.6,1.0]) #("AR1",[self.langd,0.9,0.0])#("binom",[self.langd,0.6,1.0])#("binom2",[self.langd,0.6,1.0])
        self.rankMatris=[]
        self.apparentMatris=[[j for j in i] for i in self.matris]
        self.imputationMatris=[]

    def start(self,typ, args):
        print "Into start: AR1, binom, binom2:"
        if typ=="AR1":
            fu=[ar1,args]
        elif typ=="binom":
            fu=[binom,args]
        elif typ=="binom2":
            fu=[binom2,args]
        for i in range(self.nr):
            self.matris.append(fu[0](*args))

    def rank(self): #
        kolumner=[[self.apparentMatris[i][k] for i in range(self.nr) if self.apparentMatris[i][k]!='-'] for k in range(self.langd)]
        kolumner2=[[self.apparentMatris[i][k] for i in range(self.nr)] for k in range(self.langd)]
        sortedd=[list(set([j for j in i])) for i in kolumner]#tar bort duplicates
        for i in sortedd:
            i.sort()
            
        ranger=[[] for i in range(self.nr)]
        for i in range(self.nr):
            for j in range(self.langd):
                if self.apparentMatris[i][j]=='-':
                    ranger[i].append('-')
                    continue
                k=self.apparentMatris[i][j]
                ranger[i].append(float(sortedd[j].index(k)))
        self.rankMatris=ranger
                
    def loseDataCompletelyAtRandom(self,chance):
        print "into loseDataCompletelyAtRandom "
        self.apparentMatris=[[j for j in i] for i in self.matris]
        for i in range(self.nr):
            if ran.random()<chance:
                k=ran.randint(1,self.langd-1)
                for j in range(k,self.langd):
                    self.apparentMatris[i][j]='-'

    def loseDataAtRandom(self,chance):
        print "into loseDataAtRandom "
        self.apparentMatris=[[j for j in i] for i in self.matris]
        self.rank()

        r=[max([self.rankMatris[i][j] for i in range(self.nr)]) for j in range(self.langd)]
        for i in range(len(self.matris)):
            for j in range(len(self.matris[i])):
                if self.rankMatris[i][j] not in [r[j],r[j]-1]:#[0,1]:
                    continue
                if ran.random()<chance:
                    for k in range(j+1,self.langd):
                        self.apparentMatris[i][k]='-'
                break

    def loseDataNotAtRandom(self,chance):
        print "into loseDataNotAtRandom "
        self.apparentMatris=[[j for j in i] for i in self.matris]
        self.rank()
        r1=[i[-1] for i in self.rankMatris]
        r=float(max(r1))
        for i in range(len(self.matris)):
            r1=(r-self.rankMatris[i][-1])/r
            if ran.random()<chance*r1:
                k=ran.randint(1,self.langd-1)
                for j in range(k,self.langd):
                    self.apparentMatris[i][j]='-'

    def sqDist(self):
        print "In i sqDist"
        t=[]
        langdimp=len(self.imputationMatris)
        for i in range(len(self.imputationMatris)):
            t.append((self.imputationMatris[i][-1]-self.matris[i][-1])**2)
            print "i, ",i,", sqDist: ImputationMatris[i][-1], ", self.imputationMatris[i][-1], ", - self.matris[i][-1]**2, ",  self.matris[i][-1], ",**2"
        print "sammanräknad sqDist för alla individer: ,",  float(sum(t))
        #return float((sum(t))**0.5)
        return float((sum(t)/langdimp)**0.5)
        
    def absDist(self):
        print "Into absDist"
        t=[]
        langdimp=len(self.imputationMatris)
        for i in range(len(self.imputationMatris)):
            t.append(abs(self.imputationMatris[i][-1]-self.matris[i][-1]))
            print "i ,",i,", absDist: ABS(ImputationMatris[i][-1] ,", self.imputationMatris[i][-1], ", - self.matris[i][-1], ",  self.matris[i][-1], ", )"
        print "sammanräknad absDist för alla individer: ,",  float(sum(t))
        #return float(sum(t))
        return float(sum(t)/langdimp) 

    def bias(self):
        print "Into bias"
        t=[]
        langdimp=len(self.imputationMatris)
        for i in range(len(self.imputationMatris)):
            t.append(self.imputationMatris[i][-1]-self.matris[i][-1])
            #print "i, ",i,", bias: ImputationMatris[i][-1] ,", self.imputationMatris[i][-1], ",- self.matris[i][-1], ",  self.matris[i][-1]
        #print "sammanräknad bias för alla individer: ,",  float(sum(t))
        #return float(sum(t))
        return float(sum(t)/langdimp)

    def LVCF(self):
        print "Int LVCF: "
        #Stegar igenom alla individer för att leta efter '-' på sista raden
        self.imputationMatris=[[i for i in j] for j in self.apparentMatris]
        for i in range(len(self.imputationMatris)):
            #Här trodde jag att datorn skulle sortera om så att återstående '-' kommer sist. 
            if self.imputationMatris[i][-1]=='-':
                #Jag undrar om inte tanken här är att stega igenom alla '-'. JO  SÅ ATT DEN HITTAR SENASTE OBS FÖR ATT IMPUTERA MED.
                for j in self.imputationMatris[i][::-1]:
                    if j!='-':
                        self.imputationMatris[i][-1]=j 
                        print "LVCF imputering"
                        #for p in self.imputationMatris: print "imputationMatris: ",p                        
                        break
        
            
    def lRtCF2(self, meanTrue=False, order=1):
        print "Into lRtCF2: "
        r=[]
        for i in range(len(self.rankMatris[0])):
            r.append(max([j[i] for j in self.rankMatris if j[i]!='-']))
        #print r
        impMatris=[]
       
        for i in range(len(self.rankMatris)):
            lista=[]
            for j in range(len(self.rankMatris[i])):
                if self.rankMatris[i][j]!='-':
                    rd=rateDistribution(self.rankMatris[i][j],r[j])
                    lista.append(rd)
                else:
                    break
            likelyRating=lista[-1]
            impMatris.append(likelyRating.impOrder())
            #impMatris.append(likelyRating.expectedValue())
        ''' 
        x=[]
        y=[]
        ko=[]

        for i in range(len(self.apparentMatris)):
            if self.apparentMatris[i][-1]!='-':
                ko.append([impMatris[i],self.apparentMatris[i][-1]])
        ko.sort()
        for i in ko:
            y.append(i[1])
            x.append(i[0])
            
        x=list(set(x))
        y=list(set(y))

        x.sort()
        y.sort()

        tck = InterpolatedUnivariateSpline(x, y, k=order)
        xnew = np.linspace(0,1,500)
        ynew = tck(xnew)
        '''
        self.imputationMatris=[[i for i in j] for j in self.apparentMatris]
        me=self.mean()
        meImp=self.mean(imputationTrue=True)
        for i in range(len(self.imputationMatris)):
            if self.apparentMatris[i][-1]=='-':
      ##################################          
                individList= [row[-1] for row in self.apparentMatris][::1] # printa observationer för samtliga individer vid mättillfälle z
                bortfallFinns = '-' in individList
                minOV=min(individList)
                n=len(individList)
                #Nedan: Hitta största observationen under mättillfället (max(maxOV) trots att det finns '-' i listan.
                numbers = []
                for item in individList:
                    try:
                        numbers.append(float(item))
                    except ValueError:
                            # ignore items which aren't integers
                        pass
                maxOV= max(numbers)
                #print "ind: ", i, "i: ", i," n: ", n ,"max: ", maxOV, " minOV: ", minOV
                # HAR NU GÅTT MÄTT MAXOV, MINOV, N VID ETT MÄTTILLFÄLLE
                #if self.apparentMatris[i][-1]=='-':        #Kollar bortfall för i:te personen och sista mättillfället.  
                #print "ind: ", i, "i: ", i," n: ", n ,", max: ", maxOV, ", minOV: ", minOV, ", impRank: ", likelyRating.impRank(), ", expectedValue: ", likelyRating.expectedValue()
                if bortfallFinns is True:
                    #DET FINNS MELLAN 1 OCH 9 MEN BORDE VARA 1 TILL 10
                    orderstat=(impMatris[i])/(n+1)   
                    #print "r[i]: ", r[i]                    
                    #(i+1)/float(n+1)
                    #orderstat=likelyRating.expectedValue()
                    imputering=orderstat*(maxOV-minOV)+minOV #Omvandla till ursprungligt värdeområde
                    #print self.imputationMatris[i][-1]                                          
                    self.imputationMatris[i][-1]=imputering                                    
                    #print self.imputationMatris[1][-1]
                    #print "LGCF imputering"
                    #for p in self.imputationMatris: print "imputationMatris: ",p                
      #######################################          
                '''
                if meanTrue: (DÄR KANSKE HAR ATT GÖRA MED VAR IMPUTERINGEN SKER  -I MITTEN ELLER I ÄNDARNA...)
                    #tck interpolerar så den får urspringligt värdeområde den vägen
                    #Jag måste i steg ett läggain värde omvandlade till ursprungligt värdeutrymme här
                    #steg 2 se hur jag kan få den att imputera alla mättillfällen och inte bara det sista (om jag hinner)                    
                    self.imputationMatris[i][-1]=float(tck(impMatris[i]))-meImp+me
                    print "LGCF imputering"
                    for p in self.imputationMatris: print "imputationMatris: ",p                        
                else:
                    self.imputationMatris[i][-1]=float(tck(impMatris[i]))
                    print "LGCF imputering"
                    for p in self.imputationMatris: print "imputationMatris: ",p        
                '''



    def lPCF(self, meanTrue=False, order=1):
        '''
        print "lPCF "
        r=[]
        for i in range(len(self.rankMatris[0])):
            r.append(max([j[i] for j in self.rankMatris if j[i]!='-']))
 #          r.append(max([j[i] for j in self.rankMatris]))
        #print "r: ", r
        #print "max(r): ",max(r)
        impMatris=[]
       
        for i in range(len(self.rankMatris)):
            lista=[]
            for j in range(len(self.rankMatris[i])):
                if self.rankMatris[i][j]!='-':
                    rd=self.rankMatris[i][j]/float(r[j])
                    lista.append(rd)
                else:
                    break
            print "lista: ", lista
            impMatris.append(lista[-1])
        x=[]
        y=[]
        ko=[]
        for p in self.rankMatris: print "rankMatrs: ", p
        for p in impMatris: print "impMatris: ", p
        #for p in self.apparentMatris: print "apparentMatris: ",p
        
        for i in range(len(self.apparentMatris)):
            if self.apparentMatris[i][-1]!='-':
                ko.append([impMatris[i],self.apparentMatris[i][-1]])
        ko.sort()
        for i in ko:
            y.append(i[1])
            x.append(i[0])
            
        x=list(set(x))
        y=list(set(y))

        x.sort()
        y.sort()
        tck = InterpolatedUnivariateSpline(x, y, k=order)
        xnew = np.linspace(0,1,500)
        ynew = tck(xnew)

        self.imputationMatris=[[i for i in j] for j in self.apparentMatris]
        me=self.mean()
        meImp=self.mean(imputationTrue=True)
        for i in range(len(self.imputationMatris)):
            if self.apparentMatris[i][-1]=='-':
                if meanTrue:
                    self.imputationMatris[i][-1]=float(tck(impMatris[i]))-meImp+me
                    print "LPCF imputering"
                    for p in self.imputationMatris: print "imputationMatris: ",p                    
                else:
                    #tck interpolerar så den får urspringligt värdeområde den vägen
                    #Jag måste i steg ett läggain värde omvandlade till ursprungligt värdeutrymme här (NV*(OVmax-OVmin))+OVmin
                    #steg 2 se hur jag kan få den att imputera alla mättillfällen och inte bara det sista (om jag hinner)
                    self.imputationMatris[i][-1]=float(tck(impMatris[i]))
                    print "LPCF imputering"
                    for p in self.imputationMatris: print "imputationMatris: ",p
                    '''
##############################
        r=[]
        self.imputationMatris=[[i for i in j] for j in self.apparentMatris]

        ovLast=0
        zerooneLast=0
        minovLast=0
        maxovLast=0
        imputering=0
    
            #print range(len(self.imputationMatris)) #Detta ger antalet individer!!
        for i in range(len(self.imputationMatris)): #i är individ

            if self.apparentMatris[i][-1]=='-':
                ##### Här sätta NV utifrån apparentMatris[i][j-x] 
    
                #nedan stegar sig bak längs MT från det sista tills j ger en observation istället för '-'.
                j=-1
                while self.apparentMatris[i][j]=='-':
    
                    j-=1
    
                ovLast=self.apparentMatris[i][j]
                #print obs #senaste observationen som inte var '-'
                #Här är tanken att NV baseras på ovLast den sista observationen innan bortfallen inträder. Jag vill alltså ha max och min från det MT:t.
                obsLista=[]
                obsLista=[w[j] for w in self.apparentMatris] # if i[j]!='-'
            #print "rad 444 - i: ",i
    
                numbers=[]
                for item in obsLista:
                    try:
                        numbers.append(float(item))
                    except ValueError:
                            # ignore items which aren't integers
                        pass
                maxovLast=max(numbers)
                minovLast=min(numbers)
                zerooneLast=(ovLast-minovLast)/(maxovLast-minovLast)
                leastobsLista=[]
                leastobsLista=[w[-1] for w in self.imputationMatris]
                numbers=[]
                for item in leastobsLista:
                    try:
                        numbers.append(float(item))
                    except ValueError:
                            # ignore items which aren't integers
                        pass                
                minOV=min(numbers)
                maxOV=max(numbers)             
                imputering=0 
                imputering=zerooneLast*(maxOV-minOV)+minOV

                #Lägger in imputerat värde för individ i: sista mättillfälle- Imputerade värdet beräknas utifrån  
                self.imputationMatris[i][-1]=imputering
                #for p in self.imputationMatris: print "imputationMatris: ",p

###############################
                
    def mean(self,imputationTrue=False):
        if imputationTrue:
            a=[i[-1] for i in self.imputationMatris if i[-1]!='-']
        else:
            a=[i[-1] for i in self.apparentMatris if i[-1]!='-']
        return sum(a)/float(len(a))
 
    def resultat(self): 
        print "Into resultat(self):  "
        missingness={'MCAR':self.loseDataCompletelyAtRandom,'MAR':self.loseDataAtRandom,'MNAR':self.loseDataNotAtRandom}
        measures={'bias':self.bias,'square distance':self.sqDist, 'absolute distance':self.absDist}
        imputation={'LVCF':self.LVCF,'LPCF':self.lPCF,'LGCF':self.lRtCF2}
        # interpolationsordning={1:'linear',2:'quadratic',3:'cubic'}
        interpolationsordning={1:'linear'}



        avstand=[]
        for i in missingness:
            miss=missingness[i]
        
            miss(0.25) 

            self.rank() 
            for j in ["LVCF","LPCF","LGCF"]:
                impu=imputation[j]
                if j=="LVCF":
                    #ra = 0 för LVCF
                    ra=[0]
                elif j in ['LPCF','LGCF']:
                    #ra=1,2,3 annars
                    #ra=range(1,4) 
                    ra=[1]
               
                for k in ra: 
                    if k:
                        impu(order=k)
                  
                    else:      
                        
                        impu() 
                    tempAvstand=[i,j,k]
                    tempDict={}
                    for m in ["bias","square distance","absolute distance"]:
                        tempDict[m]=measures[m]()
                    tempAvstand.append(tempDict)
                    #for p in tempAvstand: print p
                    avstand.append([t for t in tempAvstand])
        #print "De olika avstånden:"
        #for p in avstand: print p #Ger listor med det slutliga excelresultatet
        #print "avstånd: ", avstand
        return avstand


def totResultat(antal=1000):     
    print "Into totResultat():  "
    measures=["bias","square distance","absolute distance"]
    resultat={}
    r2={}

    for i in ['MCAR','MAR','MNAR']:      
#Resultat och r2.
        resultat[i]={}                   
        r2[i]={}
        for j in ['LVCF','LPCF','LGCF']:
            resultat[i][j]={}
            r2[i][j]={}

            if j=="LVCF":
                ra=[0]
            elif j in ['LPCF','LGCF']:
                #ra=range(1,4)
                ra=[1]
            for k in ra:
                resultat[i][j][k]={}
                r2[i][j][k]={}
                for m in measures:
                    resultat[i][j][k][m]=[]
                    r2[i][j][k][m]={}
                    r2[i][j][k][m]['mean']=0.0
                    r2[i][j][k][m]['standard dev']=0.0

    i=0

    while i<antal:                      

        try:                            

            l=dataStruct(personer,mattillfallen)    

            lista=l.resultat()                      
        except:
            continue
        for j in lista:
            for k in measures:
#Vad fylls resultat. resultat 
                resultat[j[0]][j[1]][j[2]][k].append(j[3][k]) 
     
        i+=1                                         


    r2={}

    for i in ['MCAR','MAR','MNAR']:             
        r2[i]={}
        for j in ['LVCF','LPCF','LGCF']:
            r2[i][j]={}
            if j=="LVCF":
                ra=[0]
            elif j in ['LPCF','LGCF']:
                #ra=range(1,4)
                ra=[1]
            for k in ra:
                r2[i][j][k]={}
                for m in measures:
                    f=sum(resultat[i][j][k][m])/float(len(resultat[i][j][k][m]))
                    sd=((sum([(ll-f)**2 for ll in resultat[i][j][k][m]])/float(len(resultat[i][j][k][m])-1.5))/len(resultat[i][j][k][m]))**0.5
                    r2[i][j][k][m]={}
                    r2[i][j][k][m]['mean']=round(f,4)
                    r2[i][j][k][m]['standard dev']=round(sd,4)
    #print len(resultat[i][j][k][m])
    #for p in r2: print p
    return r2
    
#Skriver resultatet till fil
def outputResultat(antal,namn):  
    print "Into outputResultat():  "
    measures=["bias","square distance","absolute distance"]
    dictt={'mean':'','standard dev':'std. dev.'}
    r=totResultat(antal)

#skriv over tidigare filer med samma namn
    m=open(namn + '.csv','w')
#se ovan
    m.close()                

    m=open(namn + '.csv','a')


    for i in ['MCAR','MAR','MNAR']:                    

        m.write(','+i+'\n,')                           

        for j in ['LPCF','LGCF','LVCF']:                
            if j=="LVCF":                               
                ra=[0]
            elif j in ['LPCF','LGCF']:
                #ra=range(1,4)
                ra=[1]
            for k in ra:
                if k==0:

                    m.write(j+'\n')                     
                else:
#Om LPCF el LGCF (ra=1-4) skriv + ip + k
                    #m.write(j+'ip'+str(k)+',') 
                    m.write(j+',') 

        for t in measures:                              

            for s in ['mean','standard dev']:           
#Skriv bias etc + "" om mean eller std. dev om standard deviation.
                m.write(t+' '+dictt[s]+',')             
#<
                for j in ['LPCF','LGCF','LVCF']:        
#j=LPCF etc
                    if j=="LVCF":                       
                        ra=[0]
                    elif j in ['LPCF','LGCF']:
                        #ra=range(1,4)
                        ra=[1]
                    for k in ra:                        
                        if k==0:
                            m.write(str(r[i][j][k][t][s])+'\n')  
                        else:
                            m.write(str(r[i][j][k][t][s])+',')
    #for p in r: print p 
    return r

typ='AR1'
personer=10
mattillfallen=15               
#filnamn="{}_med_{}_individer_{}_mattillfallen".format(typ,personer,mattillfallen)
#r=outputResultat(personer,filnamn) #.format(personer,mattillfallen)
#o=ar1(10,3) 
#print(o)
#print("..........................")
#i=choose(30,40)
#print(i)
#print("----------------------------")

#q=rateDistribution(50,20)
#w=q.expectedValue
#print("q = ", q)
#print("************************")
#print("w = ",w)
#x=q.expectedValue
#print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
#print(x)
#print "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV"
#w=dataStruct()
#print(w.rank)

#class rateDistribution:
#    def __init__(self,rank,antal):
#        self.rank=float(rank)
#        self.antal=float(antal)
#    def cumulativeDensity(self,y,expectedvalue=0):
#        konstant=self.rank*choose(self.antal,self.rank)
#        summa=0
#        for i in range(int(self.antal-self.rank+1)):
#            summa+=((-1)**i)*choose(self.antal-self.rank,i)*(y**(i+self.rank+expectedvalue))/float(i+self.rank+expectedvalue)
#        return konstant*summa
#    def expectedValue(self):
#        return self.rank/(self.c+1)

#g=rateDistribution()
#h=g.expectedValue(5,10)
#print(h)

#r=outputResultat(2000,'Binom_med_15_individer_10_mattillfallen')
#r=outputResultat(2000,'AR1_med_15_individer_10_mattillfallen')
r=outputResultat(10000,'10000_binom2_med_10_individer_15_mattillfallen_ver_2_output4')
