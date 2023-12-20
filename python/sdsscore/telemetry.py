"""
routines for APOGEE telemetry data
   load() reads APOGEE telemetry files and load into database
   mkhtml() reads from database and makes interactive bokeh plots
"""

import argparse
import os
import tempfile
from astropy.time import Time
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import pdb
import subprocess
from . import database
from holtztools import plots

def load(file,obs='apo',skip=30,dir='./') :
    """ read telemetry file and load to database

    Parameters
    ----------
    file : str
           Name of telemetry file to load
    obs : str, default='apo'
           apo | lco to determine sensors in file and database table
    skip : int, default=30
           load every skip entries from input file (which polls every 10 seconds)
    dir : str
          parent directory name for file

    Returns
    -------
    Table : table of mjd,quantity,value

    """

    # get rid of bad lines (NF!=6) and output to temporary file
    tmpfile = tempfile.mktemp()
    subprocess.run("awk 'NF==6' {:s}/{:s} > {:s}".format(dir,file,tmpfile), shell=True)

    # read file into table and add columns mjd, cal, quantity
    try :
        tab=Table(ascii.read(tmpfile))[::skip]
    except :
        return
    finally :
        os.remove(tmpfile)

    tab['mjd'] = Time(np.char.add(np.char.add(tab['col1'],' '),tab['col2']),format='iso').mjd
    tab['cal'] = 0.
    tab['quantity'] = '             '

    # get calibrated data for each sensor of interest
    if obs == 'apo' :
        sensors = ['E2A0820-1','TPG261-Vacuum','VacPressure', 'CCPressure',
                   "E1BEE20_1-1", "E1BEE20_1-2", "E1BEE20_0-3", "LS340-TempA",
                   "LS340-TempB", "LS340-TempC", "0_0-0", "0_1-0", "0_2-0",
                   "0_2-1", "0_2-2", "0_2-3", "0_3-2",
                   "0_0-1", "0_0-2", "0_0-3"]
    else :
        sensors = [ 'LS218ATEMP', '3424-01', '3424-02', '3424-03',
                    '3424-04', '3424-05', '3424-05', '3424-08', '3424-09',
                    'AUX3424-0', 'VacPressure', 'CCPressure', 'TPG261-Vacuum']

    for sensor in sensors :
        j = np.where(np.char.find(tab['col5'],sensor) >= 0)[0]
        print(sensor,len(j))
        quantity,cal = getcal(sensor,tab['col6'][j].astype(float),obs=obs)
        tab['quantity'][j] = quantity
        tab['cal'][j] = cal

    # load entries for sensors of interest
    tab.remove_columns(['col1','col2','col3','col4','col5','col6'])
    gd = np.where(tab['quantity'] != '             ')[0]
    d=database.DBSession()
    d.ingest('apogee_telemetry.{:s}'.format(obs),tab[gd],onconflict='update')
    d.close()

    return tab

def getcal(sensor, value, obs = 'apo') :
    """ Compute calibrated data for specified input sensor

    Parameters
    ----------
    sensor : str
             Sensor name in telemetry file
    value : array-like
             values from telemetry file
    obs : str, default='apo'
          apo | lco determines which sensors and calibration values

    Returns
    -------
    str : name corresponding to sensor
    array-like : array of calibrated values

    """

    if obs == 'apo' :
      if "E2A0820-1" in sensor :
        quantity = 'ln2'
        cal = (value-0.588)/0.02352 

      elif "TPG261-Vacuum" in sensor :
        quantity = 'vac261'
        cal = value

      elif "VacPressure" in sensor :
        quantity = 'vacMKD'
        cal = value

      elif "CCPressure" in sensor :
        quantity = 'MKS_CC'
        cal = value

      elif "E1BEE20_1-1" in sensor :
        quantity = 'Tent'
        cal = (0.9971-value)*377.87+78

      elif "E1BEE20_1-2" in sensor :
        quantity = 'DetPole'
        cal = (0.9973-value)*379.83+78

      elif "E1BEE20_0-3" in sensor :
        quantity = 'CameraAft'
        cal = (0.9941-value)*382.665+78

      elif "LS340-TempA" in sensor :
        quantity = 'DetA'
        cal = value

      elif "LS340-TempB" in sensor :
        quantity = 'DetB'
        cal = value

      elif "LS340-TempC" in sensor :
        quantity = 'DetC'
        cal = value

      elif "0_1-0" in sensor :
        quantity = 'CPMid'
        cal =  (0.9967-value)*379.14+78

      elif "0_2-0" in sensor :
        quantity = 'CPHang'
        cal = (0.9986-value)*377.87+78

      elif "0_2-1" in sensor :
        quantity = 'CPCorner'
        cal = (0.9956-value)*379.30+78

      elif "0_2-2" in sensor :
        quantity = 'Coll'
        cal = (0.9933-value)*383.50+78

      elif "0_2-3" in sensor :
        quantity = 'RadE'
        cal = (0.9963-value)*379.50+78

      elif "0_3-2" in sensor :
        quantity = 'TBoard'
        cal = (0.9963-value)*381.01+78-217 

      elif "0_0-0" in sensor :
        quantity = 'VPH'
        cal =  (0.9955-value)*380.33+78

      elif "0_0-1" in sensor :
        quantity = 'CamFwd'
        cal =  (0.9980-value)*378.33+78

      elif "0_0-2" in sensor :
        quantity = 'CamMid'
        cal = (1.0000-value)*378.20+78

      elif "0_0-3" in sensor :
        quantity = 'CamAft'
        cal =  (0.9941-value)*382.665+78


    elif obs == 'lco' :
      if "LS218ATEMP" in sensor :
        quantity = 'DetT'
        cal = value
 
      elif "3424-01" in sensor :
        quantity = 'CamFwd'
        cal=(0.4312018563-value)*453.5483725+78

      elif "3424-02" in sensor :
        quantity = 'CamMid'
        cal = (0.4311467087-value)*453.1854964+78

      elif "3424-03" in sensor :
        quantity = 'CamAft'
        cal = (0.4312375198-value)*453.7772908+78

      elif "3424-04" in sensor :
        quantity = 'CPM'
        cal = (0.4310049009-value)*453.3089586+78

      elif "3424-05" in sensor :
        quantity =  'Tent'
        cal = (0.431192618-value)*453.8104834+78

      elif "3424-05" in sensor :
        quantity = 'DPB'
        cal = (0.4310591901-value)*453.6686918+78

      elif "3424-08" in sensor :
        quantity = 'CPH'
        cal = (0.4309846833-value)*453.4413494+78

      elif "3424-09" in sensor :
        quantity = 'CPC'
        cal = (0.4309548678-value)*453.4724764+78

      elif "AUX3424-0" in sensor :
        quantity = 'ln2'
        cal = value*26+51.8 

      elif "VacPressure" in sensor :
        quantity = 'MKS'
        cal = value

      elif "CCPressure" in sensor :
        quantity = 'MKS_CC'
        cal = value

      elif "TPG261-Vacuum" in sensor :
        quantity = 'TPG'
        cal = value

    return quantity, cal

def mkhtml(obs='apo',outfile=None,skip=1) :
    """ Make interactive plots for HTML page

    Parameters
    ----------
    obs : str, apo|lco
          observatory, determines names of telemetry values
    outfile : str, default=None
          output HTML file (otherwise opens in browser)
    skip : int, default=1
          plot every skip entries
    """

    # list of tab names and quantities to plot in each tab
    if obs == 'apo' :
        names = ['ln2','Camera','ColdPlate','Vacuum']
        categories = [['ln2'],['CamFwd','CamMid','CamAft','DetA','DetB','DetC','Coll'],['CPMid','CPHang','CPCorner','VPH','Tent','RadE'],['TPG','MKS','MKS_CC']]
    else :
        names = ['ln2','Cam','Vacuum','CP','Det']
        categories = [['ln2'],['CamFwd','CamMid','CamAft'],['Vac261','MKS','MKS_CC'],['CPM','CPH','CPC'],['DetT']]

    # open database session
    d=database.DBSession()
    tabs=[]
    for category,name in zip(categories,names) :
        if name == 'Vacuum' : log=True
        else : log=False
        fig,ax = plots.bokeh_multi(1,len(category),sharex=True,sharey=True,width=1000,height=250,ylog=log,tab=name)
        for iy,quantity in enumerate(category) :
            # query database for this quantity
            out=d.query(sql="select * from apogee_telemetry.{:s} where quantity = '{:s}'".format(obs,quantity))
            if len(out) > 0 :
                plots.bokeh_plotp(ax[iy,0],out['mjd'][::skip],out['cal'][::skip],size=1,xt='MJD',yt=quantity)
            else :
                print('no data for: ',quantity)
        tabs.append(fig)

    plots.bokeh_show(tabs,tab=True,outfile=outfile)

    # close database connection
    d.close()
