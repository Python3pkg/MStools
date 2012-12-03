#!/opt/bin/python
import sys
import pymzml
from pyMSA import parseFeatureXML, getWindow
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D,art3d
from matplotlib import cm
from numpy import array
from matplotlib.collections import LineCollection
from matplotlib.colors import colorConverter,LinearSegmentedColormap
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import pyteomics.mass as mass
import subprocess
import optparse,re 
try:
    import Unimod.unimod #need to write this
except:
    print "No Unimod database module available"

def main():
    parser=optparse.OptionParser()
    parser.add_option("-m", "--mzml", dest="mzml", help="Input mzML file")
    parser.add_option("-f", "--feature", dest="feature", help="comma separated list of featureXML files")
    parser.add_option("-o", "--output", dest="outfile", help="output filename root (default='plot') For digestion pruducts the files will be called plot_1.pdf, plot_2.pdf etc.", default="plot")
    parser.add_option("-p", "--peptide", dest="peptide", help="peptide sequence")
    parser.add_option("-P", "--protein", dest="protein", action="store_true", default=False, help="the peptide sequence should be digested (a separate output file will be written for each digestion product)")
    
    parser.add_option("-e", "--enzyme", dest="enzyme", help="Enzyme to digest with.")
    parser.add_option("-c", "--modifications", dest="modifications", help="post translational modifications, comma separated as 'pos modification (aa)' ie. '3 Phospho (ST)")
    parser.add_option("-s", "--charge", dest="charge", help="charge states to consider", default='2,3,4')

    (options,args)=parser.parse_args()
    print "options settings:"
    stdval=dir(optparse.Values)
    for p in dir(options):
      if p not in stdval:
        print "--%s:\t%s"%(p,getattr(options,p))
    #print "%s %s %s"%(options, featcols, featfiles)
    #mzlist=ms1pep.listmz("DAVIDMARTIN", modifications=['9 Phospho (T)',])
    #xics=ms1pep.getXIC("JG-C1-1.mzML", mzlist)
    #ms1pep.plotXIC(xics, outfile="ms1pep.pdf", title="XIC for Peptide DAVIDMARTIN [9 Phospho (T)]") 
    peplist=[]
    enz="1"
    enzymes={ "Trypsin": "1",
             "Lys-C":"2",
             "Arg-C":"3",
             "Asp-N":"4",
             "V8-bicarb":"5",
             "V8-phosph":"6",
             "Chymotrypsin":"7",
             "CNBr":"8"
             }

    if getattr(options,"protein")==True:
        if options.enzyme !=None:
            try:
                enz=str(int(options.enzyme))
            except:
                try:
                    enz=enzymes[options.enzyme]
                except:
                    raise Exception("Invalid value for --enzyme: %s"%options.enzyme)
        peplist=proteindigest(options.peptide, enzyme=enz)
    else:
        peplist=[{"sequence": options.peptide, "start": 1, "end": len(peptide)},] 
    
    index=0   
    for p in peplist:
        index = index +1
        mods=[]
        for m in options.modifications:
            (pos, mod)=m.split(None,1)
            mp = 1+int(pos)-int(m['start'])
            if mp>=m['start'] and mp <=m['end']:
                mods.append("%s %s"%(mp,mod))
        outfile="%s.pdf"%options.outfile
        if len(peplist)>1:
            outfile="%s_%s.pdf"%(options.outfile, index)
        xics=plotXIC(getXIC(mzml, listmz(p['sequence'], charges=options.charges, modifications=mods), minrt=minrt, maxrt=maxrt), outfile=outfile, title="%s %s-%s %s"%(p['sequence'], p['start'], p['end'], ";".join(mods)))  
    

def plotXIC(xics, colours=['r','g','b','m','c','k'], outfile="plot.pdf", title="XIC plot"):
    '''Plots the given array of XIC using matplotlib. Colours can be specified as an array of codes. Default ['r','g','b','m','c','k']'''
    labels=xics.keys()
    labels.remove('rt')
    args=[]
    index=0
    maxrt=max(array(xics['rt']))
    rtdif=maxrt-xics['rt'][0]
    maxic=0
    for k in sorted(labels, key=float):
        if max(array(xics[k])) >maxic:
            maxic=max(array(xics[k]))
    for k in sorted(labels, key=float):
        args.extend([xics['rt'],array(xics[k])+(index*maxic*0.07)+(maxic*0.01),"%s-"%colours[index%len(colours)]]) 
        index=index+1
    lines=plt.plot(*args)
    plt.setp(lines, lw=0.1)
    plt.title(title)
    plt.xlabel("Retention time")
    plt.ylabel("XIC")
    index=0
    for t in sorted(labels, key=float):
        if len(t) >7:
            t=t[:7]
        plt.text(xics['rt'][0]+(0.01*rtdif), maxic*((0.07*index)+0.01), t, color=colours[index])
        index=index+1
    plt.savefig(outfile)
    plt.close()

    
def listmz(peptide, charges=[2,3,4], fixedmods={},modifications=[]):
  ''' Calculates the mz values for a given peptide with modifications for each of the charges listed in charges
    Default is to calculate 2+, 3+ and 4+
    listmz(peptide, charges=[2,3,4], fixedmods={"C": 56.0987}, modifications=['3 Phospho (STY)']
    '''
  hmass=float(Unimod.unimod.database.get_element('H')['mono_mass'])
  mz=float(mass.calculate_mass(peptide))
  for p in modifications:
      m=re.match(r'(\d+) +(.*) +\(([^\)]*)\) *$', p)
      if m:
          pos=int(m.group(1))
          label=m.group(2)
          aa=m.group(3)
          if peptide[pos-1] in aa:
              mz = mz + float(Unimod.unimod.database.get_label(label)['delta_mono_mass'])
  for k in fixedmods.keys():
      for a in peptide:
          if a==k:
              mz = mz + float(fixedmods[k])
  mzcalc=[]
  for c in charges:
      mzcalc.append(hmass+(mz/c))
  return mzcalc

def digestprotein(peptide, enzyme=1, overlap=True, unfavoured=False):
  '''Digest a protein sequence using the EMBOSS Digest program. enzyme should be one of the codes used in the digest -menu option
  return value is an array with start, end and sequence
  [{'start': xx, 'end': xx, 'sequence': 'ACDEFGHIK'},]
  '''
  digargs=["digest", "asis::%s"%peptide, "-stdout", "-menu", str(enzyme), "-auto","-rformat", "excel"]
  frags=[]
  if overlap:
      digargs.append('-overlap')
  if unfavoured:
      digargs.append('-unfavoured')
  dig=subprocess.Popen(digargs,  stdout=subprocess.PIPE)
  if dig == None:
    raise Exception("Error running digest - please check path")
  digout=dig.communicate()[0].split("\n")[1:]
  #print digout
  for d in digout:
    if len(d)>10 and d[:7] !='SeqName':
      #print d
      f=d.split()
      #print f
      frags.append({"start":f[1], "end": f[2], "sequence": peptide[int(f[1]) - 1 : int(f[2]) ]})
  return frags

def digestproteindb(peptidefile, enzyme=1, overlap=True, unfavoured=False):
  '''Digest a protein sequence using the EMBOSS Digest program. enzyme should be one of the codes used in the digest -menu option
  return value is an array with start, end and sequence
  [{'start': xx, 'end': xx, 'sequence': 'ACDEFGHIK', 'proteinID': 'P12345'},]
  '''
  digargs=["digest", peptidefile, "-stdout", "-menu", str(enzyme), "-auto","-rformat", "excel"]
  frags=[]
  if overlap:
      digargs.append('-overlap')
  if unfavoured:
      digargs.append('-unfavoured')
  dig=subprocess.Popen(digargs,  stdout=subprocess.PIPE)
  if dig == None:
    raise Exception("Error running digest - please check path")
  digout=dig.communicate()[0].split("\n")[1:]
  #print digout
  for d in digout:
    if len(d)>10 and d[:7] !='SeqName':
      #print d
      f=d.split()
      #print f
      frags.append({"start":f[1], "end": f[2], "proteinID": f[0], "sequence": peptide[int(f[1]) - 1 : int(f[2]) ]})
  return frags


    
def getXIC(mzmlfile, mzlist, tolerance=0.01, minrt=0, maxrt=None):
    '''Extract summed XIC for each  of the mz (range +/- tolerance) in mzlist optionally between minrt and maxrt (in minutes)
    returns: 
    { 'rt': [rt,...], 'mz1':[xic, ...],'mz2':[xic, ...]}

    '''   
    try:
      run=pymzml.run.Reader(mzmlfile)
    except Exception, e:
      raise Exception("Error trying to opn mzML file %s: %s"%(mzmlfile, e))
    xics={'rt':[]}
    mzkeys=[]
    for m in mzlist:
        mzkeys.append("%s"%m)
        xics["%s"%m]=[]
    if len(mzlist)==0:
      raise Exception('No mz values specified')
    scan=run.next()
    while scan['scan time'] < minrt:
      scan=run.next()
    while scan.has_key('scan time'):
      if maxrt != None and scan['scan time'] > maxrt:
        break
      if scan['ms level'] == 2:
        pass
      else:
        #print scan['scan time']
        xics['rt'].append(float(scan['scan time']))
        mzindex=0
        mzint=0
        for p in scan.peaks:
          if p[0] >mzlist[mzindex]+tolerance:
            xics[mzkeys[mzindex]].append(mzint)
            mzindex=mzindex +1
            mzint=0
          if mzindex==len(mzlist):
            break 
          if p[0] > mzlist[mzindex]-tolerance:
            mzint=mzint+p[1]
      scan=run.next()
    
    return xics
    
    
    
    
    
    
    
if __name__=="__main__":
  main()
