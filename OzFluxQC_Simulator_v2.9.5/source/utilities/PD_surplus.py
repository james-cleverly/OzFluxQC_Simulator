import sys, os

sys.path.append(os.path.abspath('scripts'))

import qcio
import numpy
import qcutils
import ast
import meteorologicalfunctions as mf

InLevel = 'L1'
cf = qcio.load_controlfile(path='controlfiles')

NEE = ast.literal_eval(cf['operations']['NEE'])
dER = ast.literal_eval(cf['operations']['dER'])
Sws = ast.literal_eval(cf['operations']['Sws'])
SR = 0.02

ds = qcio.xl_read_series(cf,InLevel)
w,f,a = qcutils.GetSeriesasMA(ds,'Uz')
c,f,a = qcutils.GetSeriesasMA(ds,'Cc')
Ah,f,a = qcutils.GetSeriesasMA(ds,'Ah')
Tv,f,a = qcutils.GetSeriesasMA(ds,'Tv')
ps,f,a = qcutils.GetSeriesasMA(ds,'ps')

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

print " Cov(wc)up: "+str(cov_up_wc)+"\n Cov(wc)down: "+str(cov_down_wc)+"\n Cor(cq)up: "+str(corr_up_cq)+"\n Cor(cq)down: "+str(corr_down_cq)+"\n"

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
print " GPP*up: "+str(GPPstar_up0)+"\n GPP*down: "+str(GPPstar_down0)+"\n GPP**up: "+str(GPPstarstar_up0)+"\n GPP**down: "+str(GPPstarstar_down0)+"\n ER*up: "+str(ERstar_up0)+"\n ER*down: "+str(ERstar_down0)+"\n ER**up: "+str(ERstarstar_up0)+"\n ER**down: "+str(ERstarstar_down0)+"\n"

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
print " alpha: "+str(alpha)+"\n net:  "+str(net)+"\n gross:  "+str(gross)+"\n"
print " alpha_hat: "+str(alphahat)+"\n ER*hat:  "+str(ERstarhat)+"\n GPP*hat:  "+str(GPPstarhat)+"\n"
print " PD case 1: "+str(PD_1)+"\n PD case 2:  "+str(PD_2)+"\n ER case 1:  "+str(ER_1)+"\n ER case 2:  "+str(ER_2)+"\n GPP case 1:  "+str(GPP_1)+"\n GPP case 2:  "+str(GPP_2)+"\n"
print filename+"\n"
print " NEE:  "+str(NEE)+"\n dER:  "+str(dER)+"\n  PD:  "+str(PD)+"\n  ER:  "+str(ER)+"\n GPP:  "+str(GPP)+"\n"

