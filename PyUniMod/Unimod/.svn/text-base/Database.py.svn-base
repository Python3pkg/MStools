import xml.etree.ElementTree as et

import sys,re

class Database():
    '''Wrapper for the Unimod database'''
    xmlns='{http://www.unimod.org/xmlns/schema/unimod_2}'
    unimodfile='unimod.xml'
    hidden=True

    def __init__(self, **kwargs):
        
        self.unimodfile=kwargs.get("file","unimod.xml")
        self.hidden=kwargs.get("hidden",True)
        node=et.parse(self.unimodfile)
        root=node.getroot()
        self.elements={}
        self.residues={}
        self.labels={}
        self.modifications={}
        self._get_elements(node)
        self._get_modifications(node)

    def _get_elements(self, node):
        for e in node.findall('%selements/%selem'%(self.xmlns,self.xmlns)):
            ea=e.attrib
            self.elements[ea['title']]=ea
            if re.match(r'[A-Z]',ea['title'][:1]):
                self.elements["%s%s"%(int(round(float(ea['mono_mass']))),ea['title'])]=ea

    def _get_modifications(self,node):
        for e in node.findall('%smodifications/%smod'%(self.xmlns, self.xmlns)):
            ma=e.attrib
            d=e.findall("%sdelta"%self.xmlns)[0]
            for k in d.attrib.keys():
                ma['delta_%s'%k]=d.attrib[k]
            ma['sites']={}
            ma['spec_group']={}
            for r in e.findall('%sspecificity'%self.xmlns):
                if self.hidden==True or r.attrib['hidden']==False:
                    ma['sites'][r.attrib['site']]=r.attrib
                    ma['sites'][r.attrib['site']]['NeutralLoss']=[]
                    # add NeutralLoss
                    for n in r.findall('%sNeutralLoss'%self.xmlns):
                        ma['sites'][r.attrib['site']]['NeutralLoss'].append(n.attrib)
                    # add to aa mods list.

                    if self.residues.has_key(r.attrib['site']):
                        self.residues[r.attrib['site']].append(ma['title'])
                    else:
                        self.residues[r.attrib['site']]=[ma['title'],]

                    if ma['spec_group'].has_key(r.attrib['spec_group']):
                        ma['spec_group'][r.attrib['spec_group']].append(r.attrib['site'])
                    else:
                        ma['spec_group'][r.attrib['spec_group']]=[r.attrib['site'],]
            self.modifications[ma['title']]=ma
                        

    def get_label(self, label):
        mod=self.modifications.get(label, None)
        return mod

    def get_element(self, name):
        el=self.elements.get(name, None)
        return el


    def list_labels(self, search):
        labels=[]
        lre=re.compile(search)
        for k in self.modifications.keys():
            l=lre.search(k)
            if l !=None:
                labels.append(k)
        return labels
    
    def get_neutral_loss(self, label, site):
        mod=self.modifications.get(label,None)
        if mod !=None:
            try:
                nl=[]
                for n in mod['sites'][site]['NeutralLoss']:
                    if n['composition']!='0':
                        nl.append(n)
                return nl
            except:
                return []
        return []

    def get_delta_mono(self, label):
        mod=self.modifications.get(label,None)
        if mod !=None:
            try:
                val=float(mod['delta_mono_mass'])
                return val
            except:
                pass
