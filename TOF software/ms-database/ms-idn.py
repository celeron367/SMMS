# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import numpy as np
import os
from scipy.signal import convolve
import epics
from epics import caget,caput,camonitor
import time


class sprctro_compair():
    thresh=3
    def __init__(self,peak_stat):
        self.peak_stat = peak_stat
        self.matched= pd.DataFrame({'m/z':[]})
        self.unmatched= pd.DataFrame({'m/z':[]})
        self.others= pd.DataFrame({'m/z':[]})
        self.matched_peak=[]
        self.score=0
    def match(self,x,y):
        if (len(self.peak_stat)==0):
            return -1
        else:
            matched = False
            for i in range(len(self.peak_stat)):            
                if (abs(self.peak_stat.iloc[i]['m/z']-x)<self.thresh):
                    #print('point matched',x,self.peak_stat.iloc[i]['m/z'],self.peak_stat.iloc[i]['count'])
                    self.matched.loc[len(self.matched)]=x
                    self.matched_peak.append(i)
                    self.score=self.score+self.peak_stat.iloc[i]['count']
                    matched=True
                    break
            if(matched==False):
                #print('no matched point')
                self.others.loc[len(self.unmatched)]=x
    def match_peaks(self,peaks):
        for i in range(len(peaks)):
            self.match(peaks.loc[i,0],peaks.loc[i,1])
        self.unmatched['m/z']=self.peak_stat.iloc[list(set([i for i in range(len(self.peak_stat))])-set(self.matched_peak))]['m/z']
        return self.score

class reference_compair():
    thresh=2
    def __init__(self,peak_stat):
        self.peak_stat = peak_stat
        self.matched= pd.DataFrame({'m/z':[]})
        self.unmatched= pd.DataFrame({'m/z':[]})
        self.others= pd.DataFrame({'m/z':[]})
        self.matched_peak=[]
        self.score=0
    def match(self,x,y):
        if (len(self.peak_stat)==0):
            return -1
        else:
            matched = False
            for i in range(len(self.peak_stat)):            
                if (abs(self.peak_stat.iloc[i]['m/z']-x)<self.thresh):
                    #print('point matched',x,self.peak_stat.iloc[i]['m/z'],self.peak_stat.iloc[i]['count'])
                    self.matched.loc[len(self.matched)]=x
                    self.matched_peak.append(i)
                    self.score=(self.peak_stat.iloc[i]['intensity']+y)+10/(abs(self.peak_stat.iloc[i]['m/z']-x)+0.01)
                    matched=True
                    break
            if(matched==False):
                #print('no matched point')
                self.others.loc[len(self.unmatched)]=x
    def match_peaks(self,peaks):
        for i in range(len(peaks)):
            self.match(peaks.loc[i,0],peaks.loc[i,1])
        self.unmatched['m/z']=self.peak_stat.iloc[list(set([i for i in range(len(self.peak_stat))])-set(self.matched_peak))]['m/z']
        return self.score

def find_peak(y):
    Y=-1*y
    kernel = [1,-1]
    dY = convolve(Y, kernel, 'valid') 

    #Checking for sign-flipping
    S = np.sign(dY)
    ddS = convolve(S, kernel, 'valid')
    #These candidates are basically all negative slope positions
    #Add one since using 'valid' shrinks the arrays
    candidates = np.where(dY < 0)[0] + (len(kernel)+1)

    #Here they are filtered on actually being the final such position in a run of
    #negative slopes
    peaks = sorted(set(candidates).intersection(np.where(ddS == 2)[0]+1 ))
    #If you need a simple filter on peak size you could use:
    alpha =-15
    peaks = np.array(peaks)[Y[peaks] < alpha]
    return peaks    
    
def get_refences(peak_df,kind,extrace):
    buff=peak_df[(peak_df['btype']==kind) & (peak_df['extrace']==extrace)]
    buff['ID']=buff['ID'].astype(np.int)
    index=set(buff['ID'])
    return buff,index


def process():
    print('开始')
    a='开始检索'
    caput('ID:INFO',a.decode('utf-8').encode('cp936'))
    global in_out_data
    global super_spectro
    global peak_df
    
    #caput('ID:STATUS',1)
    
    x=caget('MS:X')
    y=caget('MS:Y')
    
    
    ty=[]
    sc=[]
    for tp in set(super_spectro['type']):
        compair=sprctro_compair(super_spectro[super_spectro['type']==tp])
        peaks=find_peak(y)
        peaks_table=pd.DataFrame()
                #print(x[peaks])
                #print(y[peaks])
        peaks_table[0]=x[peaks]
        peaks_table[1]=y[peaks]
        score=compair.match_peaks(peaks_table)
        ty.append(tp)
        sc.append(score)
        print(tp.decode('utf-8').encode('cp936'),score)
        a='种类:'
        b='得分:'
        caput('ID:INFO',a.decode('utf-8').encode('cp936')+tp.decode('utf-8').encode('cp936')+','+b.decode('utf-8').encode('cp936')+ str(score))       
    idn_rst=ty[sc.index(max(sc))]
    print(idn_rst.decode('utf-8').encode('cp936'))
    a='鉴定结果： '
    caput('ID:INFO',a.decode('utf-8').encode('cp936')+idn_rst)

    s1=pd.Series(ty,index=sc)
    s1=s1.sort_index(ascending=False)  
    s2=pd.Series(sc,index=ty)
 #-----------------------------------------------------------------------------
    a='搜索参考谱......'
    caput('ID:INFO',a.decode('utf-8').encode('cp936'))
    idn_rst=s1.iloc[0]
    refences,index=get_refences(peak_df,idn_rst,'提取')     
    FitScore=0
    bestFit=0
    for i in index:
        ref=refences[refences['ID']==i]
        cmp=reference_compair(ref)
        ref_score=cmp.match_peaks(peaks_table)
        print(i,ref_score)
        if (ref_score>FitScore):
            bestFit=i
            FitScore=ref_score
    

    xr=refences[refences['ID']==bestFit]['m/z']
    yr=refences[refences['ID']==bestFit]['intensity']
    xc=pd.to_numeric(in_out_data.iloc[bestFit]['m/z'].replace('[','').replace(']','').split(','))
    yc=pd.to_numeric(in_out_data.iloc[bestFit]['intensity'].replace('[','').replace(']','').split(','))
        #--------------Write data to IOC----------------------------------
    caput('ID1:X',xr)
    caput('ID1:Y',yr)
    caput('ID1:REF:X',xc)
    caput('ID1:REF:Y',yc)
    caput('ID1:NAME',idn_rst.decode('utf-8').encode('cp936'))
    caput('ID1:SCORE',s2.loc[idn_rst])
#----------------------------------------------------------------------------------  
    idn_rst=s1.iloc[1]
    refences,index=get_refences(peak_df,idn_rst,'提取')     
    FitScore=0
    bestFit=0
    for i in index:
        ref=refences[refences['ID']==i]
        cmp=reference_compair(ref)
        ref_score=cmp.match_peaks(peaks_table)
        print(i,ref_score)
        if (ref_score>FitScore):
            bestFit=i
            FitScore=ref_score
    

    xr=refences[refences['ID']==bestFit]['m/z']
    yr=refences[refences['ID']==bestFit]['intensity']
    xc=pd.to_numeric(in_out_data.iloc[bestFit]['m/z'].replace('[','').replace(']','').split(','))
    yc=pd.to_numeric(in_out_data.iloc[bestFit]['intensity'].replace('[','').replace(']','').split(','))
        #--------------Write data to IOC----------------------------------
    caput('ID2:X',xr)
    caput('ID2:Y',yr)
    caput('ID2:REF:X',xc)
    caput('ID2:REF:Y',yc)
    caput('ID2:NAME',idn_rst.decode('utf-8').encode('cp936'))
    caput('ID2:SCORE',s2.loc[idn_rst])
#-------------------------------------------------------------------------------------
    idn_rst=s1.iloc[2]
    refences,index=get_refences(peak_df,idn_rst,'提取')     
    FitScore=0
    bestFit=0
    for i in index:
        ref=refences[refences['ID']==i]
        cmp=reference_compair(ref)
        ref_score=cmp.match_peaks(peaks_table)
        print(i,ref_score)
        if (ref_score>FitScore):
            bestFit=i
            FitScore=ref_score
    

    xr=refences[refences['ID']==bestFit]['m/z']
    yr=refences[refences['ID']==bestFit]['intensity']
    xc=pd.to_numeric(in_out_data.iloc[bestFit]['m/z'].replace('[','').replace(']','').split(','))
    yc=pd.to_numeric(in_out_data.iloc[bestFit]['intensity'].replace('[','').replace(']','').split(','))
        #--------------Write data to IOC----------------------------------
    caput('ID3:X',xr)
    caput('ID3:Y',yr)
    caput('ID3:REF:X',xc)
    caput('ID3:REF:Y',yc)
    caput('ID3:NAME',idn_rst.decode('utf-8').encode('cp936'))
    caput('ID3:SCORE',s2.loc[idn_rst])

    caput('ID:STATUS',0)
    a='结束'
    caput('ID:INFO',a.decode('utf-8').encode('cp936'))
    
    
    
    
def process_callback(pvname=None, value=None, host=None, **kws):
    print('process data')
    
    
    
    

global in_out_data
global super_spectro
global peak_df

if __name__=="__main__":
    global in_out_data
    global super_spectro
    global peak_df
    print('Loading database.')
    f = open('./predict_result1.csv')
    in_out_data = pd.read_csv(f,usecols=[0,1,2,3],header=None,names=['btype','extrace','m/z','intensity',])  
    in_out_data.loc[in_out_data['extrace'].str.contains('直涂'),'extrace']='直涂'
    in_out_data.loc[in_out_data['extrace'].str.contains('提取'),'extrace']='提取'   
    print('Loading index')
    f = open('./super_spectro.csv')
    super_spectro = pd.read_csv(f,usecols=[0,1,2,3,4,5,6],header=None,names=['type','extrace','m/z','count','intensity','intensity_min','intensity_max'])   
    print('Loading references data')
    f = open('./peak_list.csv')
    peak_df = pd.read_csv(f,usecols=[1,2,3,4,5],names=['ID','btype','extrace','m/z','intensity'])  
    #---------------------------监视IOC变化，等候指令----------------
    x=caget('MS:X')
    y=caget('MS:Y')
    #camonitor('MS:START',callback=process_callback)
    caput('MS:START',0)
    while(True):
        flag=caget('MS:START')
        if(flag==1):
            process()
            caput('MS:START',0)
            print('done')
        time.sleep(0.5)
        


        
    
    

    
