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
import meteorologicalfunctions as mf

log = logging.getLogger('partition.ts')

def conditional_correlation(cf,ds):
    NEE = ast.literal_eval(cf['operations']['NEE'])
    dER = ast.literal_eval(cf['operations']['dER'])
    Sws = ast.literal_eval(cf['operations']['Sws'])
    SR = 0.02
    
    w,f = putils.GetSeriesasMA(ds,'Uz')
    c,f = putils.GetSeriesasMA(ds,'Cc')
    Ah,f = putils.GetSeriesasMA(ds,'Ah')
    Tv,f = putils.GetSeriesasMA(ds,'Tv')
    ps,f = putils.GetSeriesasMA(ds,'ps')
    
    e = mf.vapourpressure(Ah,Tv)
    mr = mf.mixingratio(ps,e)
    q = mf.specifichumidity(mr) * 1000
    
    wprime = w - numpy.mean(w)
    cprime = c - numpy.mean(c)
    qprime = q - numpy.mean(q)
    wprime_sd = wprime / numpy.std(w)
    cprime_sd = cprime / numpy.std(c)
    qprime_sd = qprime / numpy.std(q)
    nRecs = len(wprime)
    
    ## separate the vertical wind into updrafts and downdrafts
    updraft_index = numpy.where(wprime>0)[0]
    downdraft_index = numpy.where(wprime<0)[0]
    
    # set up variables
    down_wprime = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    up_wprime = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    down_cprime = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    up_cprime = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    down_wprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    up_wprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    down_cprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    up_cprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    down_qprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    up_qprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    in_variables = [wprime,cprime,wprime_sd,cprime_sd,qprime_sd]
    up_variables = [up_wprime,up_cprime,up_wprime_sd,up_cprime_sd,up_qprime_sd]
    down_variables = [down_wprime,down_cprime,down_wprime_sd,down_cprime_sd,down_qprime_sd]
    
    for z in range(len(in_variables)):
        up_variables[z][updraft_index] = in_variables[z][updraft_index]
        up_variables[z] = numpy.ma.masked_values(up_variables[z], -9999.)
        down_variables[z][downdraft_index] = in_variables[z][downdraft_index]
        down_variables[z] = numpy.ma.masked_values(down_variables[z], -9999.)
    
    up_wprime = up_variables[0]
    up_cprime = up_variables[1]
    up_wprime_sd = up_variables[2]
    up_cprime_sd = up_variables[3]
    up_qprime_sd = up_variables[4]
    down_wprime = down_variables[0]
    down_cprime = down_variables[1]
    down_wprime_sd = down_variables[2]
    down_cprime_sd = down_variables[3]
    down_qprime_sd = down_variables[4]
    
    ## calcluate eddy accumulation dead band
    beta = numpy.std(w)/(numpy.mean(up_wprime_sd)-numpy.mean(down_wprime_sd))
    
    down_pass = numpy.ma.where(down_wprime_sd<-beta)[0]
    up_pass = numpy.ma.where(up_wprime_sd>beta)[0]
    
    ## generate w'/sd, c'/sd and q'/sd without data from dead band
    down_nodead_wprime = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    up_nodead_wprime = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    down_nodead_cprime = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    up_nodead_cprime = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    down_nodead_wprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    up_nodead_wprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    down_nodead_cprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    up_nodead_cprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    down_nodead_qprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    up_nodead_qprime_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) - 9999
    nodead_in_up_variables = [up_wprime,up_cprime,up_wprime_sd,up_cprime_sd,up_qprime_sd]
    nodead_in_down_variables = [down_wprime,down_cprime,down_wprime_sd,down_cprime_sd,down_qprime_sd]
    nodead_up_variables = [up_nodead_wprime,up_nodead_cprime,up_nodead_wprime_sd,up_nodead_cprime_sd,up_nodead_qprime_sd]
    nodead_down_variables = [down_nodead_wprime,down_nodead_cprime,down_nodead_wprime_sd,down_nodead_cprime_sd,down_nodead_qprime_sd]
    
    for z in range(len(in_variables)):
        nodead_up_variables[z][up_pass] = nodead_in_up_variables[z][up_pass]
        nodead_up_variables[z] = numpy.ma.masked_values(nodead_up_variables[z], -9999.)
        nodead_down_variables[z][down_pass] = nodead_in_down_variables[z][down_pass]
        nodead_down_variables[z] = numpy.ma.masked_values(nodead_down_variables[z], -9999.)
    
    up_nodead_wprime = nodead_up_variables[0]
    up_nodead_cprime = nodead_up_variables[1]
    up_nodead_wprime_sd = nodead_up_variables[2]
    up_nodead_cprime_sd = nodead_up_variables[3]
    up_nodead_qprime_sd = nodead_up_variables[4]
    down_nodead_wprime = nodead_down_variables[0]
    down_nodead_cprime = nodead_down_variables[1]
    down_nodead_wprime_sd = nodead_down_variables[2]
    down_nodead_cprime_sd = nodead_down_variables[3]
    down_nodead_qprime_sd = nodead_down_variables[4]
    
    scales = ast.literal_eval(cf['operations']['scales'])
    n_scales = len(scales) - 1
    
    mean_up_w = numpy.zeros(n_scales,dtype=numpy.float64)
    mean_up_c = numpy.zeros(n_scales,dtype=numpy.float64)
    mean_up_c_sd = numpy.zeros(n_scales,dtype=numpy.float64)
    mean_up_q_sd = numpy.zeros(n_scales,dtype=numpy.float64)
    mean_down_w = numpy.zeros(n_scales,dtype=numpy.float64)
    mean_down_c = numpy.zeros(n_scales,dtype=numpy.float64)
    mean_down_c_sd = numpy.zeros(n_scales,dtype=numpy.float64)
    mean_down_q_sd = numpy.zeros(n_scales,dtype=numpy.float64)
    
    for j in range(n_scales):
        first = scales[j]
        second = scales[j+1] - 1
        mean_up_w[j] = numpy.mean(up_nodead_wprime[first:second])
        mean_up_c[j] = numpy.mean(up_nodead_cprime[first:second])
        mean_up_c_sd[j] = numpy.mean(up_nodead_cprime_sd[first:second])
        mean_up_q_sd[j] = numpy.mean(up_nodead_qprime_sd[first:second])
        mean_down_w[j] = numpy.mean(down_nodead_wprime[first:second])
        mean_down_c[j] = numpy.mean(down_nodead_cprime[first:second])
        mean_down_c_sd[j] = numpy.mean(down_nodead_cprime_sd[first:second])
        mean_down_q_sd[j] = numpy.mean(down_nodead_qprime_sd[first:second])
    
    diff_up_w = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    diff_up_c = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    diff_up_c_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    diff_up_q_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    diff_down_w = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    diff_down_c = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    diff_down_c_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    diff_down_q_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    sqr_up_wc = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    sqr_up_cq_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    sqr_down_wc = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    sqr_down_cq_sd = numpy.ma.zeros(nRecs,dtype=numpy.float64) -9999
    cov_up_wc = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
    corr_up_cq = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
    scale_up_cq = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
    cov_down_wc = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
    corr_down_cq = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
    scale_down_cq = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
    
    for j in range(n_scales):
        first = scales[j]
        second = scales[j+1] - 1
        diff_up_w[first:second] = up_nodead_wprime[first:second] - mean_up_w[j]
        diff_up_c[first:second] = up_nodead_cprime[first:second] - mean_up_c[j]
        diff_up_w_masked =numpy.ma.masked_values(diff_up_w, -9999.)
        diff_up_c_masked =numpy.ma.masked_values(diff_up_c, -9999.)
        sqr_up_wc[first:second] = diff_up_w_masked[first:second] * diff_up_c_masked[first:second]
        sqr_up_wc_masked = numpy.ma.masked_values(sqr_up_wc, -9999.)
        remove_masks_up1 = sqr_up_wc_masked[first:second].compressed()
        n_up1 = len(remove_masks_up1)
        cov_up_wc[j] = numpy.sum(sqr_up_wc_masked[first:second])/n_up1
    
    for j in range(n_scales):
        first = scales[j]
        second = scales[j+1] - 1
        diff_down_w[first:second] = down_nodead_wprime[first:second] - mean_down_w[j]
        diff_down_c[first:second] = down_nodead_cprime[first:second] - mean_down_c[j]
        diff_down_w_masked =numpy.ma.masked_values(diff_down_w, -9999.)
        diff_down_c_masked =numpy.ma.masked_values(diff_down_c, -9999.)
        sqr_down_wc[first:second] = diff_down_w_masked[first:second] * diff_down_c_masked[first:second]
        sqr_down_wc_masked = numpy.ma.masked_values(sqr_down_wc, -9999.)
        remove_masks_down1 = sqr_down_wc_masked[first:second].compressed()
        n_down1 = len(remove_masks_down1)
        cov_down_wc[j] = numpy.sum(sqr_down_wc_masked[first:second])/n_down1
    
    for j in range(n_scales):
        first = scales[j]
        second = scales[j+1] - 1
        diff_up_c_sd[first:second] = up_nodead_cprime_sd[first:second] - mean_up_c_sd[j]
        diff_up_q_sd[first:second] = up_nodead_qprime_sd[first:second] - mean_up_q_sd[j]
        diff_up_c_sd_masked =numpy.ma.masked_values(diff_up_c_sd, -9999.)
        diff_up_q_sd_masked =numpy.ma.masked_values(diff_up_q_sd, -9999.)
        sqr_up_cq_sd[first:second] = diff_up_c_sd_masked[first:second] * diff_up_q_sd_masked[first:second]
        sqr_up_cq_sd_masked = numpy.ma.masked_values(sqr_up_cq_sd, -9999.)
        remove_masks_up1 = sqr_up_cq_sd_masked[first:second].compressed()
        n_up1 = len(remove_masks_up1)
        std_up_c = numpy.std(up_nodead_cprime_sd[first:second])
        std_up_q = numpy.std(up_nodead_qprime_sd[first:second])
        corr_up_cq[j] = numpy.sum(sqr_up_cq_sd_masked[first:second])/(n_up1*std_up_c*std_up_q)
        if Sws < 0.05:
            if (corr_up_cq[j]) > 0:
                scale_up_cq[j] = 1
            else:
                scale_up_cq[j] = (corr_up_cq[j] + 1)
        else:
            scale_up_cq[j] = (1 + corr_up_cq[j]) / 2
    
    for j in range(n_scales):
        first = scales[j]
        second = scales[j+1] - 1
        diff_down_c_sd[first:second] = down_nodead_cprime_sd[first:second] - mean_down_c_sd[j]
        diff_down_q_sd[first:second] = down_nodead_qprime_sd[first:second] - mean_down_q_sd[j]
        diff_down_c_sd_masked =numpy.ma.masked_values(diff_down_c_sd, -9999.)
        diff_down_q_sd_masked =numpy.ma.masked_values(diff_down_q_sd, -9999.)
        sqr_down_cq_sd[first:second] = diff_down_c_sd_masked[first:second] * diff_down_q_sd_masked[first:second]
        sqr_down_cq_sd_masked = numpy.ma.masked_values(sqr_down_cq_sd, -9999.)
        remove_masks_down1 = sqr_down_cq_sd_masked[first:second].compressed()
        n_down1 = len(remove_masks_down1)
        std_down_c = numpy.std(down_nodead_cprime_sd[first:second])
        std_down_q = numpy.std(down_nodead_qprime_sd[first:second])
        corr_down_cq[j] = numpy.sum(sqr_down_cq_sd_masked[first:second])/(n_down1*std_down_c*std_down_q)
        if Sws < 0.05:
            if (corr_down_cq[j]) > 0:
                scale_down_cq[j] = 1
            else:
                scale_down_cq[j] = (corr_down_cq[j] + 1)
        else:
            scale_up_cq[j] = (1 + corr_up_cq[j]) / 2
    
    log.info(' Cov(wc)up: '+str(cov_up_wc))
    log.info(' Cov(wc)down:  '+str(cov_down_wc))
    log.info(' Cor(cq)up:  '+str(corr_up_cq))
    log.info(' Cor(cq)down:  '+str(corr_down_cq))
    
    # determine ER*up and GPP**up
    ERindex = numpy.ma.where(cov_up_wc>0)[0]
    if len(ERindex) > 0:
        ERstar_up = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
        GPPstarstar_up = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
        ERstar_up[ERindex] = scale_up_cq[ERindex]*cov_up_wc[ERindex]
        GPPstarstar_up[ERindex] = (1-scale_up_cq[ERindex])*cov_up_wc[ERindex]
        ERstar_up0 = numpy.sum(ERstar_up[ERindex])
        GPPstarstar_up0 = numpy.sum(GPPstarstar_up[ERindex])
    else:
        ERstar_up0 = 0
        GPPstarstar_up0 = 0
    
    # determine ER*down and GPP**down
    ERindex2 = numpy.ma.where(cov_down_wc>0)[0]
    if len(ERindex2) > 0:
        ERstar_down = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
        GPPstarstar_down = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
        ERstar_down[ERindex2] = scale_down_cq[ERindex2]*cov_down_wc[ERindex2]
        GPPstarstar_down[ERindex2] = (1-scale_down_cq[ERindex2])*cov_down_wc[ERindex2]
        ERstar_down0 = numpy.sum(ERstar_down[ERindex2])
        GPPstarstar_down0 = numpy.sum(GPPstarstar_down[ERindex2])
    else:
        ERstar_down0 = 0
        GPPstarstar_down0 = 0
    
    # determine GPP*up and ER**up
    GPPindex = numpy.ma.where(cov_up_wc<0)[0]
    if len(GPPindex) > 0:
        GPPstar_up = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
        ERstarstar_up = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
        GPPstar_up[GPPindex] = (1-scale_up_cq[GPPindex])*-cov_up_wc[GPPindex]
        ERstarstar_up[GPPindex] = scale_up_cq[GPPindex]*-cov_up_wc[GPPindex]
        GPPstar_up0 = numpy.sum(GPPstar_up[GPPindex])
        ERstarstar_up0 = numpy.sum(ERstarstar_up[GPPindex])
    else:
        GPPstar_up0 = 0
        ERstarstar_up0 = 0
    
    # determine GPP*down and ER**down
    GPPindex2 = numpy.ma.where(cov_down_wc<0)[0]
    if len(GPPindex2) > 0:
        GPPstar_down = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
        ERstarstar_down = numpy.ma.zeros(n_scales,dtype=numpy.float64) -9999
        GPPstar_down[GPPindex2] = (1-scale_down_cq[GPPindex2])*-cov_down_wc[GPPindex2]
        ERstarstar_down[GPPindex2] = scale_down_cq[GPPindex2]*-cov_down_wc[GPPindex2]
        GPPstar_down0 = numpy.sum(GPPstar_down[GPPindex2])
        ERstarstar_down0 = numpy.sum(ERstarstar_down[GPPindex2])
    else:
        GPPstar_down0 = 0
        ERstarstar_down0 = 0
    
    ERstar = ERstar_up0 + ERstar_down0 + ERstarstar_up0 + ERstarstar_down0
    GPPstar = GPPstar_up0 + GPPstar_down0 + GPPstarstar_up0 + GPPstarstar_down0
    
    # determine the scaling parameter between C uptake and emissions:
    # (GPP*up + GPP*down + GPP**up + GPP**down) = alpha(ER*up + ER*down + ER**up + ER**down)
    alpha = ERstar / GPPstar
    log.info(' GPP*up: '+str(GPPstar_up0))
    log.info(' GPP*down:  '+str(GPPstar_down0))
    log.info(' GPP**up:  '+str(GPPstarstar_up0))
    log.info(' GPP**down:  '+str(GPPstarstar_down0))
    log.info(' ER*up: '+str(ERstar_up0))
    log.info(' ER*down:  '+str(ERstar_down0))
    log.info(' ER**up:  '+str(ERstarstar_up0))
    log.info(' ER**down:  '+str(ERstarstar_down0))
    
    net = ERstar - GPPstar
    gross = ERstar + GPPstar
    delta = NEE - net
    if delta > 0:
        ERstarhat = ERstar + delta
        GPPstarhat = GPPstar
    else:
        ERstarhat = ERstar
        GPPstarhat = GPPstar - delta
        
    grosshat = ERstarhat + GPPstarhat
    alphahat = ERstarhat / GPPstarhat
    
    # ER0: SR+PD
    if NEE > 0:
        # case 1:  GPP from 1/alpha-hat
        GPP0_1 = NEE + (NEE/alphahat)
        ER0_1 = NEE + GPP0_1
        if ER0_1 > dER:
            ER_1 = ER0_1
            GPP_1 = GPP0_1
        else:
            ER_1 = dER
            GPP_1 = ER_1 - NEE
        PD_1 = ER_1 - dER
        # case 2:  RE from alpha-hat
        ER0_2 = NEE + (NEE*alphahat)
        if ER0_2 > dER:
            ER_2 = ER0_2
        else:
            ER_2 = dER
        GPP_2 = -(NEE-ER_2)
        PD_2 = ER_2 - dER
        # sort cases by the relative values of grosshat to NEE (large case when RE+GPP is larger than NEE
        if (grosshat > (2 * NEE)):
            if ER_1 > ER_2:
                PD = PD_1
                ER = ER_1
                GPP = GPP_1
            else:
                PD = PD_2
                ER = ER_2
                GPP = GPP_2
        else:
            if ER_1 > ER_2:
                PD = PD_2
                ER = ER_2
                GPP = GPP_2
            else:
                PD = PD_1
                ER = ER_1
                GPP = GPP_1
    else:
        # case 1:  GPP from 1/alpha-hat
        GPP0_1 = -NEE + (-NEE/alphahat)
        ER0_1 = NEE + GPP0_1
        if ER0_1 > dER:
            ER_1 = ER0_1
            GPP_1 = GPP0_1
        else:
            ER_1 = dER
            GPP_1 = ER_1 - NEE
        PD_1 = ER_1 - dER
        # case 2:  RE from alpha-hat
        ER0_2 = NEE + (NEE*alphahat)
        if ER0_2 > dER:
            ER_2 = ER0_2
        else:
            ER_2 = dER
        GPP_2 = -(NEE-ER_2)
        PD_2 = ER_2 - dER
        # sort cases by the relative values of grosshat to NEE (large case when RE+GPP is larger than NEE
        if (grosshat > -(2 * NEE)):
            if ER_1 > ER_2:
                PD = PD_1
                ER = ER_1
                GPP = GPP_1
            else:
                PD = PD_2
                ER = ER_2
                GPP = GPP_2
        else:
            if ER_1 > ER_2:
                PD = PD_2
                ER = ER_2
                GPP = GPP_2
            else:
                PD = PD_1
                ER = ER_1
                GPP = GPP_1
    
    filename = str(cf['Files']['L1']['in_filename'])
    log.info(' alpha: '+str(alpha))
    log.info(' net:  '+str(net))
    log.info(' gross:  '+str(gross))
    log.info(' ER*hat:  '+str(ERstarhat))
    log.info(' GPP*hat:  '+str(GPPstarhat))
    log.info('  PD case 1:  '+str(PD_1))
    log.info('  PD case 2:  '+str(PD_2))
    log.info('  ER case 1:  '+str(ER_1))
    log.info('  ER case 2:  '+str(ER_2))
    log.info('  GPP case 1:  '+str(GPP_1))
    log.info('  GPP case 2:  '+str(GPP_2))
    log.info(str(filename))
    log.info('   NEE: '+str(NEE))
    log.info('   dER:  '+str(dER))
    log.info('   PD:  '+str(PD))
    log.info('   GPP:  '+str(GPP))
    print "\n"
    log.info(str(nRecs))
    print "\n"
    print filename+"\n"
    print " NEE:  "+str(NEE)+"\n dER:  "+str(dER)+"\n  PD:  "+str(PD)+"\n  ER:  "+str(ER)+"\n GPP:  "+str(GPP)+"\n\n"

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
    
    ldt = ds.series['DateTime']['Data']
    IncludeList = cf['Periods'].keys()
    NumDates = len(IncludeList)
    allflag = numpy.zeros(nRecs,dtype=numpy.int32) + 1
    
    for data_month in range(NumDates):
        if data_month in analysis_months:
            analysisflag = numpy.zeros(nRecs,dtype=numpy.int32) + 1
            tempkey = str(data_month)
            IncludeDateList = ast.literal_eval(cf['Periods'][tempkey])
            try:
                si = ldt.index(datetime.datetime.strptime(IncludeDateList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            
            try:
                ei = ldt.index(datetime.datetime.strptime(IncludeDateList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            
            analysisflag[si:ei] = numpy.int32(0)
            allflag[si:ei] = numpy.int32(0)
            
            zero0index = []
            zero1index = []
            if putils.cfkeycheck(cf,Base='Coefficients',ThisOne=tempkey,key='PDzero'):
                if '0' in cf['Coefficients'][tempkey]['PDzero'].keys():
                    #pdb.set_trace()
                    PDzero0 = numpy.float64(ast.literal_eval(cf['Coefficients'][tempkey]['PDzero']['0']))
                    PDzero0a = PDzero0[0]
                    PDzero0b = PDzero0[1]
                    zero0index = numpy.where((analysisflag == 0) & ((NEE > PDzero0a) & (NEE < PDzero0b)))[0]
                if '1' in cf['Coefficients'][tempkey]['PDzero'].keys():
                    PDzero1 = numpy.float64(ast.literal_eval(cf['Coefficients'][tempkey]['PDzero']['1']))
                    PDzero1a = PDzero1[0]
                    PDzero1b = PDzero1[1]
                    zero1index = numpy.where((analysisflag == 0) & ((NEE > PDzero1a) & (NEE < PDzero1b)))[0]
            PDm3_posNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDm3_posNEE'])
            PDm2_posNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDm2_posNEE'])
            PDm1_posNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDm1_posNEE'])
            PDb_posNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDb_posNEE'])
            PDm3_negNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDm3_negNEE'])
            PDm2_negNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDm2_negNEE'])
            PDm1_negNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDm1_negNEE'])
            PDb_negNEE = ast.literal_eval(cf['Coefficients'][tempkey]['PDb_negNEE'])
            
            posmonth_index = numpy.where((analysisflag == 0) & (NEE > 0))[0]
            negmonth_index = numpy.where((analysisflag == 0) & (NEE < 0))[0]
            #pdb.set_trace()
            PD[posmonth_index] = (PDm3_posNEE * (NEE[posmonth_index]*NEE[posmonth_index]*NEE[posmonth_index])) + (PDm2_posNEE * (NEE[posmonth_index]*NEE[posmonth_index])) + (PDm1_posNEE * (NEE[posmonth_index])) + PDb_posNEE
            PD[negmonth_index] = (PDm3_negNEE * (NEE[negmonth_index]*NEE[negmonth_index]*NEE[negmonth_index])) + (PDm2_negNEE * NEE[negmonth_index]*NEE[negmonth_index]) + (PDm1_negNEE * NEE[negmonth_index]) + PDb_negNEE
            
            try:
                PD[zero0index] = 0
            except:
                log.warn('no PD zero index near NEE == 0')
            
            try:
                PD[zero1index] = 0
            except:
                log.warn('no PD zero index at extreme NEE << 0')
            
            ER_day = PD + ER_dark
            GPP = ER_day - NEE
    
    Fc_day = numpy.ma.masked_where(Fsd < 1,NEE)
    nightflag = numpy.where(Fsd < 1)[0]
    Fc_day_only = numpy.ma.masked_where(allflag == 1,Fc_day)
    GPP.mask = Fc_day_only.mask
    PD.mask = Fc_day_only.mask
    ER_day.mask = Fc_day_only.mask
    
    for i in range(nRecs):
        if Fc_day.mask[i] == True:
            if allflag[i] == 1:
                PD[i] = -9999
                GPP[i] = -9999
                ER_day[i] = -9999
                GPP.mask[i] = True
                PD.mask[i] = True
                ER_day.mask[i] = True
                GPP_flag[i] = NEE_flag[i]
                PD_flag[i] = NEE_flag[i]
                ERday_flag[i] = NEE_flag[i]
            elif i in nightflag:
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
            if allflag[i] == 1:
                PD[i] = -9999
                GPP[i] = -9999
                ER_day[i] = -9999
                GPP.mask[i] = True
                PD.mask[i] = True
                ER_day.mask[i] = True
                GPP_flag[i] = NEE_flag[i]
                PD_flag[i] = NEE_flag[i]
                ERday_flag[i] = NEE_flag[i]
            elif ER_day[i] == -9999:
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
    for ThisOne in ['ER_night','ER_dark','ER_bio','ER_day','ER','PD','GPP']:
        baddataflag = numpy.where((ds.series[ThisOne]['Data'] == -9999) & (numpy.mod(ds.series[ThisOne]['Flag'],10)==0))[0]
        ds.series[ThisOne]['Flag'][baddataflag] = numpy.int32(151)
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

