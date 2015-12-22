#    pls.py
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
    Partition v0.3 24 July 2013;

    Version History:
    <<v0.0  20 July 2011>>
    <<v0.1: 21 May 2012, stand-alone gui/app set to OzFluxQCv1.8.2.b1 standards>>
    <<v0.2: 25 June 2012, stand-alone gui/app set to OzFluxQCv2.0 standards>>
    <<v0.3: 24 July 2013, stand-alone gui/app set to OzFluxQCPlusv2.5 standards>>
"""

import sys
import logging
import ast
import constants as c
import copy
import numpy
import pio
import pts
import putils
import time
import xlrd
import logging

log = logging.getLogger('partition.ls')

def l3partition(cf,ds2):
    '''Processing OzFlux_Level3 data for partitioning'''
    # make a copy of the OzFlux_Level3 data
    ds3 = copy.deepcopy(ds2)
    ds3.globalattributes['Level'] = 'L5'
    # compute daily statistics
    if putils.cfkeycheck(cf,Base='Params',ThisOne='firstMonth'):
        M1st = int(cf['Params']['firstMonth'])
    else:
        M1st = 1
    if putils.cfkeycheck(cf,Base='Params',ThisOne='secondMonth'):
        M2nd = int(cf['Params']['secondMonth'])
    else:
        M2nd = 12
    # prep nighttime ER observations
    pts.ER_nightL3(cf,ds3,M1st,M2nd,'Fc')
    pts.LightResponseCurves(ds3)
    log.info('L3 Partitioning: All done')
    return ds3
    

def l4partition(cf,ds3):
    '''Processing OzFlux_Level4 data to partition daytime ER and GPP'''
    # make a copy of the OzFlux_Level4 data
    ds4 = copy.deepcopy(ds3)
    ds4.globalattributes['Level'] = 'L4'
    if putils.cfkeycheck(cf,Base='Params',ThisOne='Fc_in'):
        Fc_in = ast.eval_literal(cf['Params']['Fc_in'])
    else:
        Fc_in = 'Fc'
    # prep nighttime ER observations
    if (putils.cfkeycheck(cf,Base='General',ThisOne='PD') and (cf['General']['PD'] == 'True')):
        Fc_in = 'NEE'
        pts.DayPD_ERGPP_TTE(cf,ds4,Fc_in)
    else:
        pts.ER_nightL4(cf,ds4,'ER_night_gapfilled')
        if 'AliceSpringsMulga' in ds4.globalattributes['site_name']:
            pts.DayERGPP_ASM(cf,ds4,Fc_in)
        elif 'TiTreeEast' in ds4.globalattributes['site_name']:
            pts.DayERGPP_TTE(cf,ds4,Fc_in)
        else:
            log.error(' Site designation undefined in ds4.globalattributes:xl_filename')
    
    log.info('L4 Partitioning: All done')
    return ds4
    
