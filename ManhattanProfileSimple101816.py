from VSG_Module import *
from sys import argv
from os import getcwd,path
try:
    F5n=argv[2]
except:
    from Tkinter import Tk
    root=Tk()
    root.attributes("-topmost", True)
    from tkFileDialog import askopenfilename
    root.withdraw()
    PosFile1=askopenfilename(title='Positional Summary File to Open',initialdir=getcwd(), filetypes=[("Text files","PositionalSummary*.tdv")])
    root.quit()
OutFileBase1=path.basename(PosFile1).split('.')[0].strip('PositionalSummary')
StoredGranularity=10000
PlottedBinsPerPoint=10
'''Columns in Positional Summary File
0. Chr
1. BinStart
2. BinLength
3. MappableSingleReads
4. FocalRepeatSingleReads
5. ChromosomalLimitedRepeatSingleReads
6. DispersedRepeatSingleReads
7. UniqueWellPositionedReadPairSpecies
8. MeanSepUniqueWellPositionedPairs
9. StdDSepUniqueWellPositionedPairs
10. UniqueShortFragReadPairSpecies_Len<=100
11. UniqueMedFragReadPairSpecies_Len<100_Len>200
12. UniqueLongFragReadPairSpecies_Len>=200
13. SingleTagmentationCandidateReads
14. UCSCLink'''

vscale1=300
vscale2=20000
hscale1=10
vset(bg=white)
ChrList1=[]
ChrD1={}  ## Keys are chromosomes, values are total length
ChrD2={}  ## Keys are chromosomes values are vertical position
BinD1={}  ## Keys are chromosomes, values are total length
BinD2={}  ## Keys are chromosomes values are vertical position
LinkD1={}
LinkD2={}
vPos1=0
MappedColumns0={3:'Uniquely Mappable Single Reads',
                4:'Focal Repeat Single Reads',
                5:'Chromosome Limited Repeat Single Reads',
                6:'Dispersed Repeat Single Reads',
                7:'Unique Well Positioned Read Pair Species',
                10:'Unique Short Frag Read Pair Species Len<=100',
                11:'Unique Med Frag Read Pair Species 100<Len<200',
                12:'Unique Long Frag Read Pair Species 100<Len<200',
                13:'Single Tagmentation Candidate Reads'}
MappedColumns1={3:'Uniquely Mappable Single Reads',
                4:'Focal Repeat Single Reads'}
dotcolor1={3:'red',4:'blue',5:'magenta',6:'cyan',7:'green',9:'magenta',10:'orange',11:'lightgreen',12:'gray'}
chrcolor1={'chrI':'red','chrII':'cyan','chrIII':'magenta','chrIV':'blue',
           'chrV':'green','chrX':'black','chrM':'orange','Repetitive':'gray',
           'OP50_Imputed':'yellow','phiX174':'brown'}
YPosChr1={'chrI':6.0,'chrII':5.0,'chrIII':4.0,'chrIV':3.0,'chrV':2.0,'chrX':1.0,'chrM':0.0}
F1=open(PosFile1,mode='rU')
Categories1={}
F1.next()
F1.next()
for L1 in F1:
    if 'OP50' in L1 or 'phiX174' in L1:
        continue
    L2=L1.strip().split()
    if not(('Unique',L2[0]) in Categories1):
        Categories1[('Unique',L2[0])]=int(L2[3])
        Categories1[('Focal_Repeat',L2[0])]=int(L2[4])
        Categories1[('Chromosomal_Repeat',L2[0])]=int(L2[5])
        Categories1[('Dispersed_Repeat',L2[0])]=int(L2[6])
        continue
    CurChr1=L2[0]
    if not(CurChr1 in ChrD1):
        ChrD1[ CurChr1 ] = 0  ## Length
        ChrD2[ CurChr1 ] = vPos1 ## Vertical Position
        vPos1+=1
    ChrD1[ CurChr1 ] +=int(L2[2]) 
    Bin1 = int(L2[1])/(PlottedBinsPerPoint*StoredGranularity)
    SegID1=(CurChr1,Bin1)
    if not SegID1 in BinD1:
        BinD1[ SegID1 ]=(ChrD2[ CurChr1 ],Bin1)
        BinD2[ SegID1 ]={i:0 for i in MappedColumns1.keys()}
        LinkD1[ SegID1 ]={i:-1 for i in MappedColumns1.keys()}
        LinkD2[ SegID1 ]={i:'http://null' for i in MappedColumns1.keys()}
    for j in MappedColumns1:
        BinD2[SegID1][j]+=int(L2[j])
        if int(L2[j])>LinkD1[ SegID1 ][j]:
            LinkD1[ SegID1 ][j]=int(L2[j])
            LinkD2[ SegID1 ][j]=L2[-1]
            
            
F1.close()
MaxValues1={j: max( [ BinD2[SegID1][j] for SegID1 in BinD2 ] ) for j in MappedColumns1.keys()}
TotalReads1=float(sum(Categories1.values()))

for g1 in BinD1:
    vPos0,hPos1=BinD1[g1]
    if not(g1[0] in YPosChr1):
        continue
    vPos1=YPosChr1[g1[0]]
    oldV1y=0.0
    for i in MappedColumns1.keys():
        v1=BinD2[g1][i]*1.0/TotalReads1  ## sum(MaxValues1.values())
        if v1==0.0:
            continue
        newV1y=oldV1y+v1
        if g1[0]=='chrM':
            col1=chrcolor1['chrM']
        else:
            col1=dotcolor1[i]
        vrect(xc=hPos1*hscale1,xr=4,
              y1=vPos1*vscale1+oldV1y*vscale2,
              y2=vPos1*vscale1+newV1y*vscale2,
              fill=col1,
              stroke='none',
              xlink=LinkD2[g1][i].replace('=HYPERLINK','').replace('"','').replace('elegans','ce11').strip('(').strip(')'))
        if 1.0*BinD2[g1][i]/TotalReads1>=0.005:
            vtext(text="{0:.1f}".format(100.0*BinD2[g1][i]/TotalReads1),x1=hPos1*hscale1+6,yc=vPos1*vscale1+(oldV1y+v1/2)*vscale2,font='DejaVu Bold 16',color=col1)
        oldV1y=newV1y
for c1 in ChrD1:
    Reads1=Categories1[('Unique',c1)]+Categories1[('Focal_Repeat',c1)]+Categories1[('Chromosomal_Repeat',c1)]
    if c1=='':
        continue
    if not(c1 in YPosChr1):
        continue
    l1=( ChrD1[ c1 ] / (PlottedBinsPerPoint*StoredGranularity) )+1
    v1=YPosChr1[ c1 ]
    vline(y1=v1*vscale1,y2=v1*vscale1,x1=-10,x2=10+l1*hscale1,stroke='black',strokewidth=4)
    vtext(text=c1+' ('+"{0:.1f}".format(100.0*Reads1/TotalReads1)+'%)',
          x1=10+l1*hscale1,
          yc=v1*vscale1,
          color='black',
          font='DejaVu Bold 48')
    for j1 in range(ChrD1[ c1 ]/1000000+1):
        vline(y1=v1*vscale1-40,y2=v1*vscale1,x1=j1*hscale1*10,x2=j1*hscale1*10,stroke='black',strokewidth=4)
        vtext(text=str(j1),y2=v1*vscale1-50,xc=j1*hscale1*10,font='DejaVu Bold 48')
        for i1 in range(1,10):
            if i1*100+j1*1000>l1:
                break
            vt1=20
            if i1==5:
                vt1=40
            vline(y1=v1*vscale1-vt1,y2=v1*vscale1,x1=j1*hscale1*10+i1*hscale1,x2=j1*hscale1*10+i1*hscale1,stroke='black',strokewidth=3)
NSomes1=['chrI','chrII','chrIII','chrIV','chrV','chrX']
X1=sum([Categories1[('Unique',y)] for y in NSomes1])
ynow1=VSG.ymin-24
vrect(x1=2,yc=ynow1,r=24,fill=dotcolor1[3],stroke=dotcolor1[3],strokewidth=1)
vtext(text='Unique Chromosomal ('+ "{0:.1f}".format(100.0*X1/TotalReads1)+'%)',x1=54,yc=ynow1,font='DejaVu Bold 48',color=black)        

X1=sum([Categories1[('Focal_Repeat',y)] for y in NSomes1])
ynow1-=60
vrect(x1=2,yc=ynow1,r=24,fill=dotcolor1[4],stroke=dotcolor1[4],strokewidth=1)
vtext(text='Local Repeat ('+ "{0:.1f}".format(100.0*X1/TotalReads1)+'%)',x1=54,yc=ynow1,font='DejaVu Bold 48',color=black)        

X1=sum([Categories1[('Chromosomal_Repeat',y)] for y in NSomes1])
ynow1-=60
vrect(x1=2,yc=ynow1,r=24,fill=dotcolor1[5],stroke=dotcolor1[5],strokewidth=1)
vtext(text='IntraChromosomal Repeat ('+ "{0:.1f}".format(100.0*X1/TotalReads1)+'%)',x1=54,yc=ynow1,font='DejaVu Bold 48',color=black)        

X1=sum([Categories1[('Dispersed_Repeat',y)] for y in NSomes1])
ynow1-=60
vrect(x1=2,yc=ynow1,r=24,fill=dotcolor1[6],stroke=dotcolor1[6],strokewidth=1)
vtext(text='Dispersed Repeat ('+ "{0:.1f}".format(100.0*X1/TotalReads1)+'%)',x1=54,yc=ynow1,font='DejaVu Bold 48',color=black)        
Types1=['Unique','Focal_Repeat','Chromosomal_Repeat','Dispersed_Repeat']

X1=sum([Categories1[(y,'chrM')] for y in Types1])
ynow1-=60
vrect(x1=2,yc=ynow1,r=24,fill=chrcolor1['chrM'],stroke=chrcolor1['chrM'],strokewidth=1)
vtext(text='Mitochondrial ('+ "{0:.1f}".format(100.0*X1/TotalReads1)+'%)',x1=54,yc=ynow1,font='DejaVu Bold 48',color=black)        



vtext(text=OutFileBase1,font='DejaVu Bold 36',color=black,x1=VSG.xmin+60,yc=VSG.ymin-60)
vdisplay('GraphicSummary_'+OutFileBase1+'_manhattan.svg')
##vclear()
##vset(bg=white)
##
##'''
##            Categories1[('Unique',L2[0])]=int(L2[4])
##            Categories1[('Focal_Repeat',L2[0])]=int(L2[5])
##            Categories1[('Chromosomal_Repeat',L2[0])]=int(L2[6])
##            Categories1[('Dispersed_Repeat',L2[0])]=int(L2[7])           
##'''
##if TotalReads1==0:
##    TotalReads1+=1
##ac1=0
##for y in ('chrI','chrII','chrIII','chrIV','chrV','chrX','chrM','OP50_Imputed','phiX174'):
##    d1=float(Categories1[('Focal_Repeat',y)])+float(Categories1[('Chromosomal_Repeat',y)])
##    if d1/TotalReads1>0.001:
##        varc(xc=0,yc=0,r=600,a0=ac1,ad=d1/TotalReads1-.001,stroke=chrcolor1[y],strokewidth=3,fill=none,
##             colorkey=y+'_FocalRepeats '+"{0:.3f}".format(d1/TotalReads1))
##    ac1+=d1/TotalReads1
##for y in ('chrI','chrII','chrIII','chrIV','chrV','chrX','chrM','OP50_Imputed','phiX174'):
##    d1=float(Categories1[('Unique',y)])
##    if d1/TotalReads1>0.001:
##        varc(xc=0,yc=0,r=600,a0=ac1,ad=d1/TotalReads1,stroke=chrcolor1[y],strokewidth=2,fill=chrcolor1[y],
##             colorkey=y+'_Unique '+"{0:.3f}".format(d1/TotalReads1))
##    ac1+=d1/TotalReads1
##varc(xc=0,yc=0,r=600,a0=ac1,ad=1.0-ac1,stroke=gray,strokewidth=2,fill=gray,colorkey="Dispersed_Repeats "+"{0:.3f}".format(1.0-ac1))
##vcolorkey()
##vdisplay()
##
