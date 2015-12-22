#    pts.py
#    Copyright (C) 2015  Dr James Cleverly
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
    Partition Data Function Module
    Used to perform the tasks queued by pls.py
    """

import sys
import ast
import constants as c
import numpy
import pio
import putils
import time
import xlrd
import xlwt
import logging
import datetime
import pdb

log = logging.getLogger('partition.ts')

def DayERGPP_ASM(cf,ds,Fc_in):
    log.info('Beginning: Daytime ER/GPP partitioning')
    Fcmg,Fc_flag = putils.GetSeriesasMA(ds,Fc_in)
    ER_umol,ER_flag = putils.GetSeriesasMA(ds,'ER_night')
    Ts,f = putils.GetSeriesasMA(ds,'Ts')
    Fsd,f = putils.GetSeriesasMA(ds,'Fsd')
    nRecs = len(Fcmg)
    
    ER_day = numpy.ma.zeros(nRecs,numpy.float64)
    
    # calculate ecosystem respiration in umol/(m2 s) modeled from light response curve (Fc-Fsd)
    Fsd_low_Ts_low_index = numpy.where((Ts < 27.5) & (Fsd < 500))[0]
    Fsd_low_Ts_high_index = numpy.where((Ts > 27.5) & (Fsd < 500))[0]
    Fsd_high_Ts_low_index = numpy.where((Ts < 36.75) & (Fsd > 500))[0]
    Fsd_high_Ts_high_index = numpy.where((Ts > 36.75) & (Fsd > 500))[0]
    
    a_low_low = 0.3687931
    b_low_low = 0.005857659
    a_low_high = 12.28563
    b_low_high = -0.06888701
    a_high_low = 0.4060834
    b_high_low = 0.06666102
    a_high_high = 29.96834
    b_high_high = -0.05043531
    
    ER_day = numpy.ma.zeros(nRecs,numpy.float64)
    ER_day[Fsd_low_Ts_low_index] = a_low_low * numpy.exp(b_low_low * Ts[Fsd_low_Ts_low_index])
    ER_day[Fsd_low_Ts_high_index] = a_low_high * numpy.exp(b_low_high * Ts[Fsd_low_Ts_high_index])
    ER_day[Fsd_high_Ts_low_index] = a_high_low * numpy.exp(b_high_low * Ts[Fsd_high_Ts_low_index])
    ER_day[Fsd_high_Ts_high_index] = a_high_high * numpy.exp(b_high_high * Ts[Fsd_high_Ts_high_index])
    
    day_index = numpy.where(Fsd > 1)
    
    Fc = ((Fcmg * (10 ** 6)) / (1000 * 44))
    Fc_day = numpy.ma.masked_where(Fsd < 1,Fc)
    GPP = numpy.ma.zeros(nRecs,numpy.float64)
    GPP.mask = Fc_day.mask
    GPP_flag = numpy.zeros(nRecs,numpy.int32)
    
    for i in range(nRecs):
        if Fc_day.mask[i] == True:
            GPP[i] = 0
            GPP.mask[i] = False
            if Fc_flag[i] == 0 or Fc_flag[i] == 10 or Fc_flag[i] == 30:
                GPP_flag[i] = numpy.int32(110)
            else:
                GPP_flag[i] = Fc_flag[i]
        else:
            if Fc_day[i] < ER_day[i]:
                GPP[i] = ER_day[i] - Fc_day[i]
                GPP.mask[i] = False
                GPP_flag[i] = 100
                ER_flag[i] = 100
            else:
                if ER_day[i] == -9999:
                    GPP[i] = -9999
                    GPP.mask[i] = True
                    GPP_flag[i] = Fc_flag[i]
                    ER_flag[i] = Fc_flag[i]
                else:
                    ER_day[i] = Fc_day[i]
                    GPP[i] = 0
                    GPP.mask[i] = False
                    GPP_flag[i] = numpy.int32(120)
                    ER_flag[i] = numpy.int32(120)
    
    ER_umol[day_index] = ER_day[day_index]
    
    putils.CreateSeries(ds,'ER_day',ER_day,Flag=0,Descr='Daytime ecosystem respiration',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'ER',ER_umol,Flag=0,Descr='Ecosystem respiration',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'GPP',GPP,Flag=0,Descr='Ecosystem Gross Primary Production',Units='umol/(m2 s)')
    
    ds.series['ER_day']['Flag'] = ER_flag
    ds.series['ER']['Flag'] = ER_flag
    ds.series['GPP']['Flag'] = GPP_flag
    
    log.info('Day ER and GPP: All done')

def DayERGPP_TTE(cf,ds,Fc_in,PD='False'):
    log.info('Beginning: ER_night and ER_dark partitioning')
    monthabr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    Fcmg,Fc_flag = putils.GetSeriesasMA(ds,Fc_in)
    Fc = Fcmg / c.Mco2
    ER_umol,ER_flag = putils.GetSeriesasMA(ds,'ER_night')
    Ts,f = putils.GetSeriesasMA(ds,'Ts')
    Fsd,f = putils.GetSeriesasMA(ds,'Fsd')
    Sws,f = putils.GetSeriesasMA(ds,'Sws')
    Hdh,f = putils.GetSeriesasMA(ds,'Hdh')
    nRecs = len(Fcmg)
    
    ER_day = numpy.ma.zeros(nRecs,numpy.float64)
    ER_dark = numpy.ma.zeros(nRecs,numpy.float64) + c.missing_value
    ER_darkFlag = numpy.ma.zeros(nRecs,numpy.float64) + 1
    
    # sort analysis into source or sink periods
    if putils.cfkeycheck(cf,Base='LRF',ThisOne='SourcePeriods'):
        ldt = ds.series['DateTime']['Data']
        IncludeList = cf['LRF']['SourcePeriods'].keys()
        NumDates = len(IncludeList)
        analysisflag = numpy.zeros(nRecs,dtype=numpy.int32) + 1
        for i in range(NumDates):
            IncludeDateList = ast.literal_eval(cf['LRF']['SourcePeriods'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(IncludeDateList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            
            try:
                ei = ldt.index(datetime.datetime.strptime(IncludeDateList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            
            analysisflag[si:ei] = numpy.int32(0)
        
        sourceindex = numpy.where(analysisflag == 0)[0]
    else:
        sourceindex = []
    
    if putils.cfkeycheck(cf,Base='LRF',ThisOne='SinkPeriods'):
        ldt = ds.series['DateTime']['Data']
        IncludeList = cf['LRF']['SinkPeriods'].keys()
        NumDates = len(IncludeList)
        analysisflag = numpy.zeros(nRecs,dtype=numpy.int32) + 1
        for i in range(NumDates):
            IncludeDateList = ast.literal_eval(cf['LRF']['SinkPeriods'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(IncludeDateList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            
            try:
                ei = ldt.index(datetime.datetime.strptime(IncludeDateList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            
            analysisflag[si:ei] = numpy.int32(0)
        
        sinkindex = numpy.where(analysisflag == 0)[0]
    else:
        sinkindex = []
    
    if len(sinkindex) != 0:
        ER_dark[sinkindex] = dark_efflux_sink(Sws[sinkindex],Ts[sinkindex],Fsd[sinkindex],Fc[sinkindex],Fc_flag[sinkindex])
        log.info('Daytime dark ER: All done for sink periods')
    
    if len(sourceindex) != 0:
        ER_dark[sourceindex] = dark_efflux_source(Sws[sourceindex],Ts[sourceindex],Fsd[sourceindex],Fc[sourceindex],Fc_flag[sourceindex])
        log.info('Daytime dark ER: All done for source periods')
    
    ER_darkFlag = numpy.ma.zeros(nRecs,numpy.float64)
    Fc_day = numpy.ma.masked_where(Fsd < 1,Fc)
    ER_dark.mask = Fc_day.mask
    for i in range(nRecs):
        if Fc_day.mask[i] == True:
            #ER_dark[i] = c.missing_value
            if Fc_flag[i] == 0 or Fc_flag[i] == 10 or Fc_flag[i] == 30:
                ER_darkFlag[i] = numpy.int32(110)
                ER_dark[i] = 0
                ER_dark.mask = False
            else:
                ER_darkFlag[i] = Fc_flag[i]
        else:
            if i in sinkindex:
                ER_darkFlag[i] = numpy.int32(130)
            elif i in sourceindex:
                ER_darkFlag[i] = numpy.int32(140)
    
    putils.CreateSeries(ds,'ER_dark',ER_dark,Flag=0,Descr='Daytime ecosystem dark respiration',Units='umol/(m2 s)')
    ds.series['ER_dark']['Flag'] = ER_darkFlag
    putils.CreateSeries(ds,'NEE',Fc,Flag=Fc_flag,Units='umol/(m2 s)')
    
    log.info('Night ER and day ER_dark: All done')

def DayPD_ERGPP_TTE(cf,ds,Fc_in):
    log.info('Beginning: Daytime PD/ER/GPP partitioning')
    
    Fsd,f = putils.GetSeriesasMA(ds,'Fsd')
    NEE,NEE_flag = putils.GetSeriesasMA(ds,Fc_in)
    ER_night,nER_flag = putils.GetSeriesasMA(ds,'ER_night')
    ER_dark,dER_flag = putils.GetSeriesasMA(ds,'ER_dark')
    Month,f = putils.GetSeriesasMA(ds,'Month')
    nRecs = len(NEE)
    MergeSeries(cf,ds,'ER_bio',[0,10,20,30,40,50,60,70,80,90,100,120,130,140,150,160,170,180,190,200])
    ER_bio,ER_bio_flag = putils.GetSeriesasMA(ds,'ER_bio')
    
    ER_day = numpy.ma.zeros(nRecs,numpy.float64) -9999
    PD = numpy.ma.zeros(nRecs,numpy.float64) -9999
    GPP = numpy.ma.zeros(nRecs,numpy.float64) -9999
    analysis_months = numpy.int32(cf['Coefficients'].keys())
    GPP_flag = numpy.zeros(nRecs,numpy.int32)
    PD_flag = numpy.zeros(nRecs,numpy.int32)
    ERday_flag = numpy.zeros(nRecs,numpy.int32)
    for data_month in [1,2,3,4,5,6,7,8,9,10,11,12,13]:
        if data_month in analysis_months:
            tempkey = str(data_month)
            PDm_negNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDm_negNEE'])
            PDb_negNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDb_negNEE'])
            PDm_posNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDm_posNEE'])
            PDb_posNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDb_posNEE'])
            GPPm_negNEE = ast.literal_eval(cf['Coefficients'][tempkey]['GPPm_negNEE'])
            GPPb_negNEE = ast.literal_eval(cf['Coefficients'][tempkey]['GPPb_negNEE'])
            GPPm_posNEE = ast.literal_eval(cf['Coefficients'][tempkey]['GPPm_posNEE'])
            GPPb_posNEE = ast.literal_eval(cf['Coefficients'][tempkey]['GPPb_posNEE'])
            GPP_zeroNEE = ast.literal_eval(cf['Coefficients'][tempkey]['GPP_zeroNEE'])
            zeroPDpos = -PDb_posNEE / PDm_posNEE
            zeroPDneg = -PDb_negNEE / PDm_negNEE
            if zeroPDpos < 0: zeroPDpos = 0
            if zeroPDneg > 0: zeroPDneg = 0
            #neg_month_day_index = numpy.where((Month != data_month) or (Fsd < 1) or (NEE > 0))[0]
            #pos_month_day_index = numpy.where((Month != data_month) or (Fsd < 1) and (NEE < 0))[0]
            zeromonth_index = numpy.where(Month == data_month)[0]
            PD[zeromonth_index] = 0
            ER_day[zeromonth_index] = ER_bio[zeromonth_index]
            GPP[zeromonth_index] = GPP_zeroNEE
            posmonth_index = numpy.where((Month == data_month) & (NEE > zeroPDpos))[0]
            PD[posmonth_index] = (PDm_posNEE * NEE[posmonth_index]) + PDb_posNEE
            ER_day[posmonth_index] = PD[posmonth_index] + ER_dark[posmonth_index]
            GPP[posmonth_index] = (GPPm_posNEE * NEE[posmonth_index]) + GPPb_posNEE
            negmonth_index = numpy.where((Month == data_month) & (NEE < zeroPDneg))[0]
            PD[negmonth_index] = (PDm_negNEE * NEE[negmonth_index]) + PDb_negNEE
            ER_day[negmonth_index] = PD[negmonth_index] + ER_dark[negmonth_index]
            GPP[negmonth_index] = (GPPm_negNEE * NEE[negmonth_index]) + GPPb_negNEE
            Fsd_index = numpy.where(Fsd < 1)
            PD[Fsd_index] = 0
            ER_day[Fsd_index] = 0
            GPP[Fsd_index] = 0
            GPP_flag[Fsd_index] = numpy.int32(110)
            PD_flag[Fsd_index] = numpy.int32(110)
            ERday_flag[Fsd_index] = numpy.int32(110)
    
    Fc_day = numpy.ma.masked_where(Fsd < 1,NEE)
    nightflag = numpy.where(Fsd < 1)[0]
    GPP.mask = Fc_day.mask
    PD.mask = Fc_day.mask
    ER_day.mask = Fc_day.mask
    for i in range(nRecs):
        if Fc_day.mask[i] == True:
            if i in nightflag:
                PD[i] = 0
                GPP[i] = 0
                ER_day[i] = 0
                GPP.mask[i] = False
                PD.mask[i] = False
                ER_day.mask[i] = False
                GPP_flag[i] = numpy.int32(110)
                PD_flag[i] = numpy.int32(110)
                ERday_flag[i] = numpy.int32(110)
            else:
                PD[i] = -9999
                GPP[i] = -9999
                ER_day[i] = -9999
                GPP.mask[i] = True
                PD.mask[i] = True
                ER_day.mask[i] = True
                GPP_flag[i] = NEE_flag[i]
                PD_flag[i] = NEE_flag[i]
                ERday_flag[i] = NEE_flag[i]
        else:
            if ER_day[i] == -9999:
                GPP[i] = -9999
                PD[i] = -9999
                GPP.mask[i] = True
                PD.mask[i] = True
                ER_day.mask[i] = True
                GPP_flag[i] = numpy.int32(151)
                PD_flag[i] = numpy.int32(151)
                ERday_flag[i] = numpy.int32(151)
            else:
                GPP.mask[i] = False
                PD.mask[i] = False
                ER_day.mask[i] = False
                GPP_flag[i] = numpy.int32(150)
                PD_flag[i] = numpy.int32(150)
                ERday_flag[i] = numpy.int32(150)
    
    putils.CreateSeries(ds,'ER_day',ER_day,Flag=0,Descr='Daytime ecosystem respiration',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'GPP',GPP,Flag=0,Descr='Ecosystem Gross Primary Production',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'PD',PD,Flag=0,Descr='Ecosystem decomposition by way of photodegradation',Units='umol/(m2 s)')
    ds.series['ER_day']['Flag'] = ERday_flag
    ds.series['PD']['Flag'] = PD_flag
    ds.series['GPP']['Flag'] = GPP_flag
    MergeSeries(cf,ds,'ER',[0,10,20,30,40,50,60,70,80,90,100,120,130,140,150,160,170,180,190,200])
    log.info('Day PD, ER and GPP: All done')

def dark_efflux_sink(Sws,Ts,Fsd,Fc,Fc_flag):
    nRecs = len(Sws)
    # calculate dark efflux in umol/(m2 s) modeled from light response curves (Fc-Fsd)
    Sws_low_Ts_low_index = numpy.where((Ts < 30) & (Sws < 0.05))[0]
    Sws_low_Ts_high_index = numpy.where((Ts > 30) & (Sws < 0.05))[0]
    Sws_high_Ts_low_index = numpy.where((Ts < 30) & (Sws > 0.05))[0]
    Sws_high_Ts_high_index = numpy.where((Ts > 30) & (Sws > 0.05))[0]
    
    a_low_low = 0.07938007
    b_low_low = 0.1168854
    a_low_high = 11.49373
    b_low_high = -0.04912147
    a_high_low = 0.6322741
    b_high_low = 0.05652428
    a_high_high = 29.28972
    b_high_high = -0.07353866
    
    ER_dark = numpy.ma.zeros(nRecs,numpy.float64)
    ER_umol = numpy.ma.zeros(nRecs,numpy.float64) + c.missing_value
    ER_dark[Sws_low_Ts_low_index] = a_low_low * numpy.exp(b_low_low * Ts[Sws_low_Ts_low_index])
    ER_dark[Sws_low_Ts_high_index] = a_low_high * numpy.exp(b_low_high * Ts[Sws_low_Ts_high_index])
    ER_dark[Sws_high_Ts_low_index] = a_high_low * numpy.exp(b_high_low * Ts[Sws_high_Ts_low_index])
    ER_dark[Sws_high_Ts_high_index] = a_high_high * numpy.exp(b_high_high * Ts[Sws_high_Ts_high_index])
    
    #pdb.set_trace()
    
    day_index = numpy.where(Fsd > 1)
    ER_umol[day_index] = ER_dark[day_index]
    return ER_umol

def dark_efflux_source(Sws,Ts,Fsd,Fc,Fc_flag):
    nRecs = len(Sws)
    a = 13.55803
    b = -0.1011087
    
    ER_dark = numpy.ma.zeros(nRecs,numpy.float64)
    ER_umol = numpy.ma.zeros(nRecs,numpy.float64) + c.missing_value
    zeroindex = numpy.where(Sws<0.00745747722936)
    ER_dark = a * Sws + b
    ER_dark[zeroindex] = 0
    
    day_index = numpy.where(Fsd > 1)
    ER_umol[day_index] = ER_dark[day_index]
    return ER_umol

def do_attributes(cf,ds):
    """
        Import attriubes in xl2nc control file to netCDF dataset.  Included
        global and variable attributes.  Also attach flag definitions to global
        meta-data for reference.
        
        Usage pts.do_attributes(cf,ds)
        cf: control file
        ds: data structure
        """
    log.info(' Getting the attributes given in control file')
    if 'Global' in cf.keys():
        for gattr in cf['Global'].keys():
            ds.globalattributes[gattr] = cf['Global'][gattr]
        ds.globalattributes['Flag000'] = 'Good data'
        ds.globalattributes['Flag001'] = 'QA/QC: -9999 in level 1 dataset'
        ds.globalattributes['Flag002'] = 'QA/QC: L2 Range Check'
        ds.globalattributes['Flag003'] = 'QA/QC: L2 Diurnal SD Check'
        ds.globalattributes['Flag004'] = 'QA/QC: CSAT Diagnostic'
        ds.globalattributes['Flag005'] = 'QA/QC: LI7500 Diagnostic'
        ds.globalattributes['Flag006'] = 'QA/QC: Excluded Dates'
        ds.globalattributes['Flag007'] = 'QA/QC: Excluded Hours'
        ds.globalattributes['Flag008'] = 'albedo: bad Fsd < threshold (290 W/m2 default) only if bad time flag not set'
        ds.globalattributes['Flag009'] = 'albedo: bad time flag (not midday 10.00 to 14.00)'
        ds.globalattributes['Flag010'] = 'Corrections: Apply Linear'
        ds.globalattributes['Flag011'] = 'Corrections/Combinations: Coordinate Rotation (Ux, Uy, Uz, UxT, UyT, UzT, UxA, UyA, UzA, UxC, UyC, UzC, UxUz, UxUx, UxUy, UyUz, UxUy, UyUy)'
        ds.globalattributes['Flag012'] = 'Corrections/Combinations: Massman Frequency Attenuation Correction (Coord Rotation, Tv_CSAT, Ah_HMP, ps)'
        ds.globalattributes['Flag013'] = 'Corrections/Combinations: Virtual to Actual Fh (Coord Rotation, Massman, Ta_HMP)'
        ds.globalattributes['Flag014'] = 'Corrections/Combinations: WPL correction for flux effects on density measurements (Coord Rotation, Massman, Fhv to Fh, Cc_7500_Av)'
        ds.globalattributes['Flag015'] = 'Corrections/Combinations: Ta from Tv'
        ds.globalattributes['Flag016'] = 'Corrections/Combinations: L3 Range Check'
        ds.globalattributes['Flag017'] = 'Corrections/Combinations: L3 Diurnal SD Check'
        ds.globalattributes['Flag018'] = 'Corrections/Combinations: u* filter'
        ds.globalattributes['Flag019'] = 'Corrections/Combinations: Gap coordination'
        ds.globalattributes['Flag021'] = 'Corrections/Combinations: Data-flag mismatch'
        ds.globalattributes['Flag022'] = 'Corrections/Combinations: L=0 in zeta=z/L'
        ds.globalattributes['Flag030'] = 'GapFilling: Gap Filled by ANN (SOLO)'
        ds.globalattributes['Flag031'] = 'GapFilling: Gap not filled'
        ds.globalattributes['Flag040'] = 'GapFilling: Gap Filled from one-minute average'
        ds.globalattributes['Flag050'] = 'GapFilling: Gap Filled from alternate'
        ds.globalattributes['Flag060'] = 'GapFilling: Gap Filled by interpolation'
        ds.globalattributes['Flag070'] = 'GapFilling: Gap Filled by replacement from paired tower'
        ds.globalattributes['Flag080'] = 'GapFilling: u* from Fh'
        ds.globalattributes['Flag081'] = 'GapFilling: u* not from Fh'
        ds.globalattributes['Flag082'] = 'GapFilling: L4 Range Check'
        ds.globalattributes['Flag083'] = 'GapFilling: L4 Diurnal SD Check'
        ds.globalattributes['Flag084'] = 'GapFilling: L5 Range Check'
        ds.globalattributes['Flag085'] = 'GapFilling: L5 Diurnal SD Check'
        ds.globalattributes['Flag086'] = 'GapFilling: L6 Range Check'
        ds.globalattributes['Flag087'] = 'GapFilling: L6 Diurnal SD Check'
        ds.globalattributes['Flag090'] = 'Partitioning Night: ER computed from exponential temperature response curves'
        ds.globalattributes['Flag100'] = 'Partitioning Day: GPP/ER computed from light-response curves, GPP = ER - Fc'
        ds.globalattributes['Flag110'] = 'Partitioning Day: GPP night mask'
        ds.globalattributes['Flag120'] = 'Partitioning Day: Fc > ER, GPP = 0, ER = Fc'
        ds.globalattributes['Flag130'] = 'Partitioning Day: ER_dark from sink period (+NEP) light response curves'
        ds.globalattributes['Flag140'] = 'Partitioning Day: ER_dark from source  period (+NEE) light response curves'
        ds.globalattributes['Flag150'] = 'Partitioning Day: PD, ER_day & GPP from conditional correlation'
        ds.globalattributes['Flag151'] = 'Partitioning Day: no solution from conditional correlation'
        ds.globalattributes['Flag161'] = 'Footprint: Date filter'
        ds.globalattributes['Flag162'] = 'Footprint: no solution'
        ds.globalattributes['Flag171'] = 'Penman-Monteith: bad rav or rSm only if bad Uavg, bad Fe and bad Fsd flags not set'
        ds.globalattributes['Flag172'] = 'Penman-Monteith: bad Fe < threshold (0 W/m2 default) only if bad Fsd flag not set'
        ds.globalattributes['Flag173'] = 'Penman-Monteith: bad Fsd < threshold (10 W/m2 default)'
        ds.globalattributes['Flag174'] = 'Penman-Monteith: Uavg == 0 (undefined aerodynamic resistance under calm conditions) only if bad Fe and bad Fsd flags not set'
        ds.globalattributes['Flag180'] = 'Penman-Monteith 2-layer: rav_base short-circuit'
        ds.globalattributes['Flag190'] = 'Penman-Monteith 2-layer: rav_top short-circuit'
        ds.globalattributes['Flag191'] = 'Penman-Monteith 2-layer: rav_top not short-circuit (rav_base undefined)'
        ds.globalattributes['Flag200'] = 'Penman-Monteith 2-layer: parallel circuit'
        ds.globalattributes['Flag201'] = 'Penman-Monteith 2-layer: not parallel circuit (rav_full short-circuit)'
        ds.globalattributes['Flag211'] = 'Bulk Richardson number flags: delta_U=0 (RiB infinite)'
    for ThisOne in ds.series.keys():
        if ThisOne in cf['Variables']:
            if 'Attr' in cf['Variables'][ThisOne].keys():
                ds.series[ThisOne]['Attr'] = {}
                for attr in cf['Variables'][ThisOne]['Attr'].keys():
                    ds.series[ThisOne]['Attr'][attr] = cf['Variables'][ThisOne]['Attr'][attr]

def do_functions(cf,ds):
    log.info(' Resolving functions given in control file')
    for ThisOne in cf['Variables'].keys():
        if 'Function' in cf['Variables'][ThisOne].keys():
            ds.series[ThisOne] = {}
            FunctionList = cf['Variables'][ThisOne]['Function'].keys()
            if len(FunctionList) == 1:
                i = 0
                if 'Square' in cf['Variables'][ThisOne]['Function'][str(i)].keys() and 'Parent' in cf['Variables'][ThisOne]['Function'][str(i)]['Square'].keys():
                    Parent = cf['Variables'][ThisOne]['Function'][str(i)]['Square']['Parent']
                    ds.series[ThisOne]['Data'] = pts.Square(ds.series[Parent]['Data'])
                    nRecs = numpy.size(ds.series[ThisOne]['Data'])
                    if 'Flag' not in ds.series[ThisOne].keys():
                        ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,numpy.int32)
                        if 'Flag' in ds.series[Parent]:
                            ds.series[ThisOne]['Flag'] = ds.series[Parent]['Flag']
                        else:
                            ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,numpy.int32)
                elif 'SquareRoot' in cf['Variables'][ThisOne]['Function'][str(i)].keys() and 'Parent' in cf['Variables'][ThisOne]['Function'][str(i)]['SquareRoot'].keys():
                    Parent = cf['Variables'][ThisOne]['Function'][str(i)]['SquareRoot']['Parent']
                    ds.series[ThisOne]['Data'] = pts.SquareRoot(ds.series[Parent]['Data'])
                    nRecs = numpy.size(ds.series[ThisOne]['Data'])
                    if 'Flag' not in ds.series[ThisOne].keys():
                        ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,numpy.int32)
                        if 'Flag' in ds.series[Parent]:
                            ds.series[ThisOne]['Flag'] = ds.series[Parent]['Flag']
                        else:
                            ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,numpy.int32)
                else:
                    log.error ('Function missing or unknown for variable'+ThisOne)
                    return
            else:
                for i in range(len(FunctionList)):
                    if 'Square' in cf['Variables'][ThisOne]['Function'][str(i)].keys() and 'Parent' in cf['Variables'][ThisOne]['Function'][str(i)]['Square'].keys():
                        Parent = cf['Variables'][ThisOne]['Function'][str(i)]['Square']['Parent']
                        ds.series[ThisOne]['Data'] = pts.Square(ds.series[Parent]['Data'])
                        nRecs = numpy.size(ds.series[ThisOne]['Data'])
                        if 'Flag' not in ds.series[ThisOne].keys():
                            ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,numpy.int32)
                            if 'Flag' in ds.series[Parent]:
                                ds.series[ThisOne]['Flag'] = ds.series[Parent]['Flag']
                            else:
                                ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,numpy.int32)
                    elif 'SquareRoot' in cf['Variables'][ThisOne]['Function'][str(i)].keys() and 'Parent' in cf['Variables'][ThisOne]['Function'][str(i)]['SquareRoot'].keys():
                        Parent = cf['Variables'][ThisOne]['Function'][str(i)]['SquareRoot']['Parent']
                        ds.series[ThisOne]['Data'] = pts.SquareRoot(ds.series[Parent]['Data'])
                        nRecs = numpy.size(ds.series[ThisOne]['Data'])
                        if 'Flag' not in ds.series[ThisOne].keys():
                            ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,numpy.int32)
                            if 'Flag' in ds.series[Parent]:
                                ds.series[ThisOne]['Flag'] = ds.series[Parent]['Flag']
                            else:
                                ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,numpy.int32)
                    else:
                        log.error ('Function missing or unknown for variable'+ThisOne)
                        return

def get_nightsums(Data):
    """
        Get nightly sums and averages on nights when no 30-min observations are missing.
        Nights with missing observations return a value of -9999
        Values returned are sample size (Num), sums (Sum) and average (Av)
        
        Usage qcts.get_nightsums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(Data.mask == False)[0]
    Num = numpy.size(li)
    if Num == 0:
        Sum = -9999
        Av = -9999
    else:
        x = 0
        for i in range(len(Data)):
            if Data.mask[i] == True:
                x = x + 1
        
        if x == 0:
            Sum = numpy.ma.sum(Data[li])
            Av = numpy.ma.mean(Data[li])
        else:
            Sum = -9999
            Av = -9999
    
    return Num, Sum, Av

def get_yearmonthdayhourminutesecond(cf,ds):
    """
        Gets year, month, day, hour, and if available seconds, from
        excel-formatted Timestamp
        
        Usage qcts.get_yearmonthdayhourminutesecond(cf,ds)
        cf: control file
        ds: data structure
        """
    log.info(' Getting date and time variables')
    # set the date mode for PC or MAC versions of Excel dates
    datemode = 0
    if cf['General']['Platform'] == 'Mac': datemode = 1
    nRecs = len(ds.series['xlDateTime']['Data'])
    Year = numpy.array([-9999]*nRecs,numpy.int32)
    Month = numpy.array([-9999]*nRecs,numpy.int32)
    Day = numpy.array([-9999]*nRecs,numpy.int32)
    Hour = numpy.array([-9999]*nRecs,numpy.int32)
    Minute = numpy.array([-9999]*nRecs,numpy.int32)
    Second = numpy.array([-9999]*nRecs,numpy.int32)
    Hdh = numpy.array([-9999]*nRecs,numpy.float64)
    Ddd = numpy.array([-9999]*nRecs,numpy.float64)
    flag = numpy.zeros(nRecs)
    for i in range(nRecs):
        DateTuple = xlrd.xldate_as_tuple(ds.series['xlDateTime']['Data'][i],datemode)
        Year[i] = numpy.float64(DateTuple[0])
        Month[i] = numpy.float64(DateTuple[1])
        Day[i] = numpy.float64(DateTuple[2])
        Hour[i] = numpy.float64(DateTuple[3])
        Minute[i] = numpy.float64(DateTuple[4])
        Second[i] = numpy.float64(DateTuple[5])
        Hdh[i] = numpy.float64(DateTuple[3])+numpy.float64(DateTuple[4])/60.
        Ddd[i] = ds.series['xlDateTime']['Data'][i] - xlrd.xldate.xldate_from_date_tuple((Year[i],1,1),datemode)
    putils.CreateSeries(ds,'Year',Year,Flag=flag,Descr='Year',Units='none')
    putils.CreateSeries(ds,'Month',Month,Flag=flag,Descr='Month',Units='none')
    putils.CreateSeries(ds,'Day',Day,Flag=flag,Descr='Day',Units='none')
    putils.CreateSeries(ds,'Hour',Hour,Flag=flag,Descr='Hour',Units='none')
    putils.CreateSeries(ds,'Minute',Minute,Flag=flag,Descr='Minute',Units='none')
    putils.CreateSeries(ds,'Second',Second,Flag=flag,Descr='Second',Units='none')
    putils.CreateSeries(ds,'Hdh',Hdh,Flag=flag,Descr='Decimal hour of the day',Units='none')
    putils.CreateSeries(ds,'Ddd',Ddd,Flag=flag,Descr='Decimal day of the year',Units='none')

def TaLightResponseCurves(ds):
    '''Prepare data into temperature classes for light response curve evaluation
    '''
    log.info('Preparing data for light response curves')
    Fcmg,f = putils.GetSeriesasMA(ds,'Fc')
    Ta,f = putils.GetSeriesasMA(ds,'Ta_HMP')
    Sws,f = putils.GetSeriesasMA(ds,'Sws')
    Fsd,f = putils.GetSeriesasMA(ds,'Fsd')
    nRecs = len(Fcmg)
    
    Fc0 = (Fcmg * (10 ** 6)) / (1000 * 44)
    Fc1 = numpy.ma.masked_where((Fsd < 1),Fc0)
    Fc_5to11 = numpy.ma.masked_where(Ta > 11,Fc1)
    Fc_11to16 = numpy.ma.masked_where((Ta < 11) | (Ta > 16),Fc1)
    Fc_16to21 = numpy.ma.masked_where((Ta < 16) | (Ta > 21),Fc1)
    Fc_21to26 = numpy.ma.masked_where((Ta < 21) | (Ta > 26),Fc1)
    Fc_26to31 = numpy.ma.masked_where((Ta < 26) | (Ta > 31),Fc1)
    Fc_31to36 = numpy.ma.masked_where((Ta < 31) | (Ta > 36),Fc1)
    Fc_36to42 = numpy.ma.masked_where((Ta < 36),Fc1)
    Fc_swc0to5 = numpy.ma.masked_where(Sws > 0.05,Fc1)
    Fc_swc5to10 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.10),Fc1)
    Fc_swc10to20 = numpy.ma.masked_where((Sws < 0.10) | (Sws > 0.20),Fc1)
    Fc_swc20to40 = numpy.ma.masked_where(Sws < 0.20,Fc1)
    
    putils.CreateSeries(ds,'Fc_5to11',Fc_5to11,FList=['Fc'],Descr='Daytime carbon flux, 5-11 C air temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_11to16',Fc_11to16,FList=['Fc'],Descr='Daytime carbon flux, 11-16 C air temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_16to21',Fc_16to21,FList=['Fc'],Descr='Daytime carbon flux, 16-21 C air temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_21to26',Fc_21to26,FList=['Fc'],Descr='Daytime carbon flux, 21-26 C air temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_26to31',Fc_26to31,FList=['Fc'],Descr='Daytime carbon flux, 26-31 C air temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_31to36',Fc_31to36,FList=['Fc'],Descr='Daytime carbon flux, 31-36 C air temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_36to42',Fc_36to42,FList=['Fc'],Descr='Daytime carbon flux, 36-42 C air temperature',Units='umol/(m2 s)')    
    putils.CreateSeries(ds,'Fc_swc0to5',Fc_swc0to5,FList=['Fc'],Descr='Daytime carbon flux, 0 to 5% soil water content',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_swc5to10',Fc_swc5to10,FList=['Fc'],Descr='Daytime carbon flux, 5 to 10% soil water content',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_swc10to20',Fc_swc10to20,FList=['Fc'],Descr='Daytime carbon flux, 10 to 20% soil water content',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_swc20to40',Fc_swc20to40,FList=['Fc'],Descr='Daytime carbon flux, 20 to 40% soil water content',Units='umol/(m2 s)')
    
    log.info('Light Response Curves: All ready')

def LightResponseCurves(ds):
    '''Prepare data into temperature classes for light response curve evaluation
    '''
    log.info('Preparing data for light response curves')
    Fcmg,f = putils.GetSeriesasMA(ds,'Fc')
    Ts,f = putils.GetSeriesasMA(ds,'Ts')
    Sws,f = putils.GetSeriesasMA(ds,'Sws')
    Fsd,f = putils.GetSeriesasMA(ds,'Fsd')
    nRecs = len(Fcmg)
    
    Fc0 = (Fcmg * (10 ** 6)) / (1000 * 44)
    Fc1 = numpy.ma.masked_where((Fsd < 1),Fc0)
    Fc_6to11 = numpy.ma.masked_where(Ts > 11,Fc1)
    Fc_11to16 = numpy.ma.masked_where((Ts < 11) | (Ts > 16),Fc1)
    Fc_16to21 = numpy.ma.masked_where((Ts < 16) | (Ts > 21),Fc1)
    Fc_21to26 = numpy.ma.masked_where((Ts < 21) | (Ts > 26),Fc1)
    Fc_26to31 = numpy.ma.masked_where((Ts < 26) | (Ts > 31),Fc1)
    Fc_31to36 = numpy.ma.masked_where((Ts < 31) | (Ts > 36),Fc1)
    Fc_36to41 = numpy.ma.masked_where((Ts < 36) | (Ts > 41),Fc1)
    Fc_41to46 = numpy.ma.masked_where((Ts < 41) | (Ts > 46),Fc1)
    Fc_46to51 = numpy.ma.masked_where((Ts < 46) | (Ts > 51),Fc1)
    Fc_51to56 = numpy.ma.masked_where((Ts < 51) | (Ts > 56),Fc1)
    Fc_56to61 = numpy.ma.masked_where((Ts < 56),Fc1)
    Fc_6to16 = numpy.ma.masked_where(Ts > 16,Fc1)
    Fc_16to26 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Fc1)
    Fc_26to36 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Fc1)
    Fc_36to46 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Fc1)
    Fc_46to56 = numpy.ma.masked_where((Ts < 46) | (Ts > 56),Fc1)
    Fc_56to66 = numpy.ma.masked_where((Ts < 56),Fc1)
    Fc_46to59 = numpy.ma.masked_where((Ts < 46),Fc1)
    Fc_swc0to5 = numpy.ma.masked_where(Sws > 0.05,Fc1)
    Fc_swc5to10 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.10),Fc1)
    Fc_swc10to20 = numpy.ma.masked_where((Sws < 0.10) | (Sws > 0.20),Fc1)
    Fc_swc20to40 = numpy.ma.masked_where(Sws < 0.20,Fc1)
    
    putils.CreateSeries(ds,'Fc_Ts6to11',Fc_6to11,FList=['Fc'],Descr='Daytime carbon flux, 6-11 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts11to16',Fc_11to16,FList=['Fc'],Descr='Daytime carbon flux, 11-16 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts16to21',Fc_16to21,FList=['Fc'],Descr='Daytime carbon flux, 16-21 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts21to26',Fc_21to26,FList=['Fc'],Descr='Daytime carbon flux, 21-26 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts26to31',Fc_26to31,FList=['Fc'],Descr='Daytime carbon flux, 26-31 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts31to36',Fc_31to36,FList=['Fc'],Descr='Daytime carbon flux, 31-36 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts36to41',Fc_36to41,FList=['Fc'],Descr='Daytime carbon flux, 36-41 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts41to46',Fc_41to46,FList=['Fc'],Descr='Daytime carbon flux, 41-46 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts46to51',Fc_46to51,FList=['Fc'],Descr='Daytime carbon flux, 46-51 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts51to56',Fc_51to56,FList=['Fc'],Descr='Daytime carbon flux, 51-56 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts56to61',Fc_56to61,FList=['Fc'],Descr='Daytime carbon flux, 56-61 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts46to59',Fc_46to59,FList=['Fc'],Descr='Daytime carbon flux, 46-59 C soil temperature',Units='umol/(m2 s)')
    
    putils.CreateSeries(ds,'Fc_Ts6to16',Fc_6to16,FList=['Fc'],Descr='Daytime carbon flux, 6-16 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts16to26',Fc_16to26,FList=['Fc'],Descr='Daytime carbon flux, 16-26 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts26to36',Fc_26to36,FList=['Fc'],Descr='Daytime carbon flux, 26-36 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts36to46',Fc_36to46,FList=['Fc'],Descr='Daytime carbon flux, 36-46 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts46to56',Fc_46to56,FList=['Fc'],Descr='Daytime carbon flux, 46-56 C soil temperature',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_Ts56to66',Fc_56to66,FList=['Fc'],Descr='Daytime carbon flux, 56-66 C soil temperature',Units='umol/(m2 s)')
    
    putils.CreateSeries(ds,'Fc_swc0to5',Fc_swc0to5,FList=['Fc'],Descr='Daytime carbon flux, 0 to 5% soil water content',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_swc5to10',Fc_swc5to10,FList=['Fc'],Descr='Daytime carbon flux, 5 to 10% soil water content',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_swc10to20',Fc_swc10to20,FList=['Fc'],Descr='Daytime carbon flux, 10 to 20% soil water content',Units='umol/(m2 s)')
    putils.CreateSeries(ds,'Fc_swc20to40',Fc_swc20to40,FList=['Fc'],Descr='Daytime carbon flux, 20 to 40% soil water content',Units='umol/(m2 s)')
    
    log.info('Light Response Curves: All ready')

def ER_nightL3(cf,ds,M1st,M2nd,Fc_in):
    log.info('Beginning: Night-time sums for Q10 response curves')
    monthabr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    Fcmg,f = putils.GetSeriesasMA(ds,Fc_in)
    Ts,f = putils.GetSeriesasMA(ds,'Ts')
    Sws,f = putils.GetSeriesasMA(ds,'Sws')
    Fsd,f = putils.GetSeriesasMA(ds,'Fsd')
    Hdh,f = putils.GetSeriesasMA(ds,'Hdh')
    Day,f = putils.GetSeriesasMA(ds,'Day')
    Month,f = putils.GetSeriesasMA(ds,'Month')
    nRecs = len(Fcmg)
    Night = numpy.ma.zeros(nRecs)
    
    Fc = ((Fcmg * (10 ** 6)) / (1000 * 44)) * 1800 / 1000
    Fc_night = numpy.ma.masked_where(Fsd > 1,Fc)
    Sws_night = numpy.ma.masked_where(Fsd > 1,Sws)
    Ts_k = Ts + c.C2K
    Ts_k_night = numpy.ma.masked_where(Fsd > 1,Ts_k)
    
    putils.CreateSeries(ds,'Fc_night',Fc_night,FList=['Fc'],Descr='Nighttime carbon flux',Units='mmol/(m2 30-min)')
    putils.CreateSeries(ds,'Sws_night',Sws_night,FList=['Sws'],Descr='Soil water content',Units='fraction')
    putils.CreateSeries(ds,'Ts_K_night',Ts_k_night,FList=['Ts'],Descr='Soil temperature',Units='K')
    
    xlFileName = cf['Files']['L3']['xlFilePath']+cf['Files']['L3']['xlFileName']
    
    xlFile = xlwt.Workbook()
    
    for i in range(nRecs):
        if Hdh[i] < 12:
            if Day[i] > 1:
                Night[i] = Day[i] - 1
            else:
                if Month[i] > 1:
                    prevMonth = Month[i] - 1
                else:
                    prevMonth = 12
                if prevMonth == 1 or prevMonth == 3 or prevMonth == 5 or prevMonth == 7 or prevMonth == 8 or prevMonth == 10 or prevMonth == 12:
                    dRan = 31
                if prevMonth == 2:
                    if ds.series['Year']['Data'][0] % 4 == 0:
                        dRan = 29
                    else:
                        dRan = 28
                if prevMonth == 4 or prevMonth == 6 or prevMonth == 9 or prevMonth == 11:
                    dRan = 30
                Night[i] = dRan
                Month[i] = prevMonth
        else:
            Night[i] = Day[i]
    
    NightList = ['Fc_night','Sws_night','Ts_K_night']
    
    for ThisOne in NightList:
        log.info(' Doing nightly sums for '+ThisOne)
        Units = ds.series[ThisOne]['Attr']['units']
        xlSheet = xlFile.add_sheet(ThisOne)
        xlRow = 0
        xlSheet.write(xlRow,4,Units)
        xlRow = xlRow + 1
        xlSheet.write(xlRow,0,'Month')
        xlSheet.write(xlRow,1,'Day')
        xlSheet.write(xlRow,2,'n')
        xlSheet.write(xlRow,3,ThisOne+'_Tot')
        xlSheet.write(xlRow,4,ThisOne+'_Av')
        data = numpy.ma.masked_where(abs(ds.series[ThisOne]['Data']-numpy.float64(-9999))<c.eps,ds.series[ThisOne]['Data'])
        
        for m in range(M1st,M2nd+1):
            if m == 1 or m == 3 or m == 5 or m == 7 or m == 8 or m == 10 or m == 12:
                dR = 31
            if m == 2:
                if ds.series['Year']['Data'][0] % 4 == 0:
                    dR = 29
                else:
                    dR = 28
            if m == 4 or m == 6 or m == 9 or m == 11:
                dR = 30
            
            for n in range(1,dR+1):
                xlRow = xlRow + 1
                ni = numpy.where((Month==m) & (Night==n) & (Fsd<1))[0]
                Num,Sum,Av = get_nightsums(data[ni])
                xlSheet.write(xlRow,0,monthabr[m-1])
                xlSheet.write(xlRow,1,n)
                xlSheet.write(xlRow,2,Num)
                xlSheet.write(xlRow,3,Sum)
                xlSheet.write(xlRow,4,Av)
    
    log.info(' Saving Excel file '+xlFileName)
    xlFile.save(xlFileName)
    
    log.info('Night sums: All done')

def ER_nightL4(cf,ds,ER_in):
    '''Ingests gap-filled nocturnal ER from seasonal/swc/temperature exponential models.'''
    log.info('Beginning: Integrating nocturnal ER')
    InLevel = 'L4ER'
    OutLevel = 'L4ER'
    pio.autoxl2nc(cf,InLevel,OutLevel)
    ds4day = pio.nc_read_series(cf,'L4ER')
    
    ER,fd = putils.GetSeriesasMA(ds4day,ER_in)
    ER_day,fd = putils.GetSeriesasMA(ds4day,'Day')
    Month_ER,fd = putils.GetSeriesasMA(ds4day,'Month')
    Fsd,f = putils.GetSeriesasMA(ds,'Fsd')
    Hdh,f = putils.GetSeriesasMA(ds,'Hdh')
    Day,f = putils.GetSeriesasMA(ds,'Day')
    Month,f = putils.GetSeriesasMA(ds,'Month')
    nRecs = len(Fsd)
    nNights = len(ER)
    ER_night = numpy.ma.zeros(nRecs,numpy.float64)
    Night = numpy.ma.zeros(nRecs)
    night_flag = numpy.ma.zeros(nRecs,numpy.int32)
    
    log.info(' Night ER: shifting night/times to center on midnight')
    for i in range(nRecs):
        #if Month[i] == 1 or Month[i] == 3 or Month[i] == 5 or Month[i] == 7 or Month[i] == 8 or Month[i] == 10 or Month[i] == 12:
        #    dRan = 31
        #if Month[i] == 2:
        #    if ds.series['Year']['Data'][0] % 4 == 0:
        #        dRan = 29
        #    else:
        #        dRan = 28
        #if Month[i] == 4 or Month[i] == 6 or Month[i] == 9 or Month[i] == 11:
        #    dRan = 30
        #    
        #if Hdh[i] < 12:
        #    if Day[i] > 1:
        #        Night[i] = Day[i] - 1
        #    else:
        #        if Month[i] > 1:
        #            prevMonth = Month[i] - 1
        #        else:
        #            prevMonth = 12
        #        Night[i] = dRan
        #        Month[i] = prevMonth
        #else:
        #    Night[i] = Day[i]
        if Hdh[i] < 12:
            if Day[i] > 1:
                Night[i] = Day[i] - 1
            else:
                if Month[i] > 1:
                    prevMonth = Month[i] - 1
                else:
                    prevMonth = 12
                if prevMonth == 1 or prevMonth == 3 or prevMonth == 5 or prevMonth == 7 or prevMonth == 8 or prevMonth == 10 or prevMonth == 12:
                    dRan = 31
                if prevMonth == 2:
                    if ds.series['Year']['Data'][0] % 4 == 0:
                        dRan = 29
                    else:
                        dRan = 28
                if prevMonth == 4 or prevMonth == 6 or prevMonth == 9 or prevMonth == 11:
                    dRan = 30
                Night[i] = dRan
                Month[i] = prevMonth
        else:
            Night[i] = Day[i]
    
    log.info(' Night ER: filling L4 night ER')
    for i in range(nRecs):
        for z in range(nNights):
            if Night[i] == ER_day[z]:
                if Month[i] == Month_ER[z]:
                    ER_night[i] = ER[z]
                    night_flag[i] = ds4day.series[ER_in]['Flag'][z]
    
    noindex = numpy.where(ER_night == 0)
    ER_night[noindex] = numpy.float64(-9999)
    night_flag[noindex] = 1
    index = numpy.where(Fsd > 1)
    ER_night[index] = numpy.float64(0)
    night_flag[index] = 110
    
    putils.CreateSeries(ds,'ER_night',ER_night,Flag=0,Descr='Nighttime ecosystem respiration',Units='umol/(m2 s)')
    ds.series['ER_night']['Flag'] = night_flag
    
    log.info('Night ER: All done')
def MergeSeries(cf,ds,series,okflags):
    """
        Merge two series of data to produce one series containing the best data from both.
        Calling syntax is: MergeSeries(cf,ds,series,okflags)
         where ds is the data structure containing all series
               series (str) is the label of the destination series
               okflags (list) is a list of QC flag values for which the data is considered acceptable
        If the QC flag for Primary is in okflags, the value from Primary is placed in destination.
        If the QC flag for Primary is not in okflags but the QC flag for Secondary is, the value
        from Secondary is placed in Destination.
        """
    # check to see if the series is specified in the control file
    section = putils.get_cfsection(cf,series=series)
    if len(section)==0: return
    # check to see if the entry for series in the control file has the MergeSeries key
    if 'MergeSeries' not in cf[section][series].keys(): return
    # now get the source list and the standard name
    srclist, standardname = putils.GetMergeSeriesKeys(cf,series,section=section)
    nSeries = len(srclist)
    if nSeries==0:
        log.info(' MergeSeries: no input series specified for '+str(series))
        return
    if nSeries==1:
        log.info(' Merging series '+str(srclist)+' into '+series)
        if srclist[0] not in ds.series.keys():
            log.error('  MergeSeries: primary input series'+srclist[0]+'not found for'+str(series))
            return
        data = ds.series[srclist[0]]['Data'].copy()
        flag = ds.series[srclist[0]]['Flag'].copy()
        attr = ds.series[srclist[0]]['Attr'].copy()
        SeriesNameString = srclist[0]
    else:
        log.info(' Merging series '+str(srclist)+' into '+series)
        if srclist[0] not in ds.series.keys():
            log.error('  MergeSeries: primary input series'+srclist[0]+'not found')
            return
        data = ds.series[srclist[0]]['Data'].copy()
        flag = ds.series[srclist[0]]['Flag'].copy()
        attr = ds.series[srclist[0]]['Attr'].copy()
        SeriesNameString = srclist[0]
        srclist.remove(srclist[0])
        for ThisOne in srclist:
            if ThisOne in ds.series.keys():
                SeriesNameString = SeriesNameString+', '+ThisOne
                indx1 = numpy.zeros(numpy.size(data),dtype=numpy.int32)
                indx2 = numpy.zeros(numpy.size(data),dtype=numpy.int32)
                for okflag in okflags:
                    index = numpy.where((flag==okflag))[0]                             # index of acceptable primary values
                    indx1[index] = 1                                                   # set primary index to 1 when primary good
                    index = numpy.where((ds.series[ThisOne]['Flag']==okflag))[0]       # same process for secondary
                    indx2[index] = 1
                index = numpy.where((indx1!=1)&(indx2==1))[0]           # index where primary bad but secondary good
                data[index] = ds.series[ThisOne]['Data'][index]         # replace bad primary with good secondary
                #pdb.set_trace()
                otherflag = numpy.int32(ds.series[ThisOne]['Flag'].copy())
                flag[index] = otherflag[index]
            else:
                log.error('  MergeSeries: secondary input series'+ThisOne+'not found')
    attr["long_name"] = attr["long_name"]+", merged from " + SeriesNameString
    putils.CreateSeries(ds,series,data,Flag=flag)

