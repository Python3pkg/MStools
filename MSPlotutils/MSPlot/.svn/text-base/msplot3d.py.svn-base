#!/opt/bin/python
import sys
try:
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
    import optparse 
except:
    raise Exception("msplot3d requires the pyMS, matplotlib and numpy modules.")

def main():
    parser=optparse.OptionParser()
    parser.add_option("-m", "--mzml", dest="mzml", help="Input mzML file")
    parser.add_option("-f", "--feature", dest="feature", help="comma separated list of featureXML files")
    parser.add_option("-o", "--output", dest="outfile", help="output filename (default='plot.pdf')", default="plot.pdf")
    parser.add_option("-s", "--showX11", dest="show", help="Show X11 interactive plot.", default=True, action="store_true")
    parser.add_option("-x", "--xrot", dest="xrot", help="Xaxis rotation for plot (default (30)", default=30, type="int")
    parser.add_option("-y", "--yrot", dest="yrot", help="Yaxis rotation for plot (default: -45)", default=-45, type="int")
    parser.add_option("-l", "--limit", dest="clip", help="Clip threshold for plot", default=500, type="int")
    parser.add_option("-c", "--colours", dest="colours", help="list of colours with which to plot features. Colours will be recycled when needed (default='r,g,m,c,y,k')", default="r,g,m,c,y,k")
    parser.add_option("-r", "--rtmin", dest="rtmin", help="minimum retention time", type="float")
    parser.add_option("-R", "--rtmax", dest="rtmax", help="maximum retention time", type="float")
    parser.add_option("-t", "--mzmin", dest="mzmin", help="minimum m/z value", type="float")
    parser.add_option("-T", "--mzmax", dest="mzmax", help="maximum m/z value", type="float")
    parser.add_option("-w", "--rtwindow", dest="rtwindow", help="Window around plot area in which to identify feature peaks", type="float", default=20.0)
    parser.add_option("-p","--plotms2", dest="plotms2",help="Plot MS2 events", action="store_true", default=True)
    parser.add_option("-d","--ms2win", dest="ms2win",help="MS2 ion capture window size", type="float", default=2.0)

    (options,args)=parser.parse_args()
    print "options settings:"
    stdval=dir(optparse.Values)
    for p in dir(options):
      if p not in stdval:
        print "--%s:\t%s"%(p,getattr(options,p))
    #  print "--%s"%p

    run=None
    featfiles=[]
    featcols=[]
    minrt=0
    maxrt=0
    minmz=0
    maxmz=0
    thresh=10
    try:
      featcols=options.colours.split(',')
      featfiles=options.feature.split(",")
    
    except:
      print "no features or no colours specified"
      featfiles=[]
      
    try:
      minrt =options.rtmin
      maxrt=options.rtmax
      minmz=options.mzmin
      maxmz=options.mzmax
      thresh=options.clip
      ms2win=options.ms2win
    except:
      raise Exception( "not all options for rt and m/z have been set")
 
    plot3d(options.mzml, featfiles=featfiles, featcols=featcols,outfile=options.outfile, show=False, xrot=options.xrot, yrot=options.yrot, minrt =options.rtmin,
           maxrt=options.rtmax,
           minmz=options.mzmin,
           maxmz=options.mzmax,
           thresh=options.clip,
           ms2win=options.ms2win,
           rtwindow=options.rtwindow,
           plotms2=options.plotms2
)

def plot3d (mzml, featfiles=[], outfile='plot.pdf', show=False, xrot=30, yrot=-45, featcols=['r','g','b','m','c','k'], thresh=100, minrt=None,maxrt=None,minmz=None, maxmz=None, ms2win=2.0, rtwindow=20.0, plotms2=True):
    try:
      run=pymzml.run.Reader(mzml, MS1_Precision = 20e-6 )
    except Exception, e:
      raise Exception( "Error opening mzML file: %s"%e)

    isline=True
    maxint=0
    mzdat=[]
    intdat=[]
    rtdat=[]
    verts=[]
    msmsruns=[]
    print featcols
    if 1:
      feat=run.next()
      while feat['scan time']*60 < minrt:
        feat=run.next()
      print "starting scan processing\n"
      while feat['scan time']*60 < maxrt:
        while feat['ms level'] != 1:
          if feat['MS:1000744'] >minmz-ms2win and feat['MS:1000744']<maxmz:  
            msmsruns.append({"total ion current": feat['total ion current'], 'parent m/z': feat['MS:1000744'], 'charge': feat['MS:1000511'],'scan time': feat['scan time']*60}) 
          feat=run.next()
        mzrow=[]
        introw=[]
        rtrow=[]
        rt=feat['scan time']*60
        #print "processing scan at %s\n"%rt
        if minrt>rt:
          minrt=rt
        if maxrt<rt:
          maxrt=rt
        #print "peaks: %s"%len(feat.peaks)
        for p in feat.peaks:
          if p[1]>thresh:
            if isline==False:
            #  #print "peak start"
              pass
            isline=True
          else:
            isline=False
            #print "peak end"
            if len(mzrow):
              #print "peak end"
              rtdat.append(rtrow)
              mzdat.append(mzrow)
              intdat.append(introw)
              mzrow=[]
              introw=[]
              rtrow=[]

          if isline and p[0]> minmz and p[0]<maxmz:
              #print "adding data"
              mzrow.append(p[0])
              introw.append(p[1])
              rtrow.append(rt)
              if maxint <p[1]:
                maxint=p[1]
          #rtrow.append(rt)
        if len(mzrow):
          mzdat.append(mzrow)
          intdat.append(introw)
          rtdat.append(rtrow)
      #rtdat.append(array(rtrow))
      #  verts.append(array(zip(mzrow, introw)))
      # print "added spectrum at rt %s\n"%rt
        feat=run.next()
    featuredata=[]
    for ff in range(len(featfiles)):
      zval=(ff+1)*10000
      print "getting features from %s"%featfiles[ff]
      features=parseFeatureXML.Reader(featfiles[ff])
      fl=getWindow.FeatureLocation(features)
      featureplot=[]
      print "Window %s,%s,%s,%s"%(minrt, maxrt, minmz, maxmz)
      for f in fl.getFeatures_rtWindow(minrt-rtwindow, maxrt+rtwindow):
        print features['mz']
        if float(features['mz'])<maxmz and float(features['mz'])>minmz:
          #print "found feature %s"%features['convexhull']
          ymin=float(features['convexhull'][0][0]['y'])
          ymax=ymin
          xmin=float(features['convexhull'][0][1]['x'])
          xmax=xmin
          for ch in features['convexhull']:
            #print "%s %s %s"%(xmax, ymax, ch)
            if float(ch[0]['y'])>ymax:
              ymax=float(ch[0]['y'])
            if float(ch[0]['y'])<ymin:
              ymin=float(ch[0]['y'])
            if float(ch[1]['x'])>xmax:
              xmax=float(ch[1]['x'])
            if float(ch[1]['x'])<xmin:
              xmin=float(ch[1]['x'])
            if ymin<minmz:
              ymin=minmz
            if xmin<minrt:
              xmin=minrt
            if ymax>maxmz:
              ymax=maxmz
            if xmax>maxrt:
              xmax=maxrt
          if xmin< maxrt and xmax >minrt:
            xs=[ymin, ymax, ymax, ymin,ymin]
            zs=[zval,zval,zval,zval,zval]
            ys=[xmin,xmin,xmax,xmax,xmin]
          #print "%s %s %s"%(xs, ys, zs)
            featuredata.append((xs,ys,zs,ff%len(featcols)))
    #poly = LineCollection(verts)
    #print "lengths: rt %s, int %s, mz %s"%(len(rtdat), len(intdat),len(mzdat))
    print "plotting figure\n"            
    plot3d_data((mzdat,rtdat,intdat), msmsruns, featuredata, outfile=outfile, show=show, xrot=xrot, yrot=yrot, featcols=featcols, thresh=thresh, ms2win=ms2win, rtwindow=rtwindow, plotms2=plotms2, bounds=(minmz,maxmz, minrt, maxrt))

def plot3d_data ( msdata, msmsruns=[], featuredata=[], outfile='plot.pdf', show=False, xrot=30, yrot=-45, featcols=['r','g','b','m','c','k'], thresh=100, ms2win=2.0, rtwindow=20.0, plotms2=True, bounds=(None, None,None,None), maxint=0, title=None):
    '''
    msdata is a tuple of arrays in the order mz, rt, intensity. msmsdata is an array of ms2 features. 
    bounds is the bounding box as a tuple (minmz, maxmz, minrt, maxrt)
    maxint is the peak intensity
    '''
    (mzdat,rtdat,intdat)=msdata
    (minmz,maxmz,minrt, maxrt) =bounds
    plt.tick_params(labelsize=8)
    if title !=None:
        plt.title(title)
    fig=plt.figure()
    ax=fig.gca(projection="3d")
    cm=LinearSegmentedColormap.from_list("scale", ((0.0,'#ccccff'),(0.5,'#0000ff'),(1.0,'k') ))
    ax.set_ylabel("RT (secs)")
    ax.set_xlabel("m/z")
    ax.set_zlabel("Intensity")
    

    if minmz !=None and maxmz !=None:
        ax.set_xlim3d(minmz, maxmz)
    if maxint >0:
        ax.set_zlim3d(10,maxint)
    if minrt !=None and maxrt !=None:
        ax.set_ylim3d(minrt,maxrt)
    ax.view_init(xrot, yrot)
    #print "%s MS2 runs seen in region: %s"%(len(msmsruns), msmsruns)
    msmsruns.reverse()

    for f in featuredata:

    #surf=ax.(array(rtdat),array(mzdat), array(intdat))
    #fig.colorbar(surf,shrink=0.5, aspect=5)
    #plt.ylim(0.1,maxint)
        print 'plotting feature at %s, %s, %s: colour %s'%(f[0],f[1],f[2],featcols[f[3]])
        po=ax.plot(f[0],f[1],f[2], featcols[f[3]])
        
    ms2index=0
    
    for d in range(len(rtdat)-1,-1, -1):
      if len(rtdat[d])==0:
          continue
      #print "RT at %i is %s"%(d, rtdat[d]) 
      relrt=((rtdat[d][0]-minrt)/(0.1+maxrt-minrt))
      if relrt<=0:
          relrt=0.0
      rg=(1.0-relrt)*0.5 * float(int(2.0*(1.0-relrt)))
      if rg<=0.0:
          rg=0.0
      #if rg>1:
      #    rg=1
      
      b=1.0 - (((2.0*relrt)-1.0)*float(int(2.0*relrt)))
      if b<=0.0:
          b=0.0
      #if b>1:
      #    b=1
      print "%s rg %s %s %s %s %s"%(len(rtdat), rg,b,rtdat[d][0],minrt, maxrt )
      try:
          #print "testing for MS2 spectra"
          while ms2index<len(msmsruns) and msmsruns[ms2index]['scan time']>rtdat[d][0]:
              # print "%s %s"%(ms2index, msmsruns[ms2index])
              mz=msmsruns[ms2index]['parent m/z']
              rt=msmsruns[ms2index]['scan time']
              tic=msmsruns[ms2index]['total ion current']
              mzl=mz-(0.2*ms2win)
              mzu=mz+(0.8*ms2win)
              if mzl < minmz:
                  mzl=minmz
              if mzu > maxmz:
                  mzu=maxmz
              if maxint >0:
                  tic=min(tic, maxint)
              xs=[mzl, mzu,mzu,mzl,mzl]
              ys=[rt,rt,rt,rt,rt]
              zs=[0,0,tic,tic,0]
              charge=int(msmsruns[ms2index]['charge'])-1
              ms2index = ms2index + 1
              if plotms2==True:
                  ax.plot(xs,ys,zs, featcols[charge])
                  #print "plotting MS2 with colour %s %s"%(charge, featcols[charge])
      except Exception,e:
          print "Error reading MS2 runs %s: %s"%(ms2index, e)
#print "rt: %s rtrel %s rg %s b %s\n"%(rtdat[d][0], relrt,rg,b)
      try:
          ints=intdat[d]
          if maxint>0:
              ints=[min(v, maxint) for v in intdat[d]]
          ax.plot(mzdat[d], rtdat[d], ints, c=(rg,rg,b), zdir='z', linewidth=0.2)
      except Exception, e:
        raise Exception("color at %s, scale %s, rg %s b %s :%s\n"%(rtdat[d][0],relrt, rg, b, e)) 

    #        feat=mpatches.Rectangle((xmin,ymin), xmax-xmin, ymax-ymin, alpha=0.5,c='r')
    #        featureplot1.append(feat)
    #collection=PatchCollection(featureplot1)
    #ax.add_collection(collection)
    #plt.savefig('defaults.pdf', format='pdf')
    for tick in ax.xaxis.get_ticklabels():
        tick.set_fontsize(10)
        tick.set_rotation('vertical')
    ax.xaxis.set_label('\n\nm/z')
    if outfile !=None:
      plt.savefig(outfile)
#    if options.show:
#      plt.show()
    #print "%s %s %s"%(options, featcols, featfiles)
if __name__=="__main__":
  main()

