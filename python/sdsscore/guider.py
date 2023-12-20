import glob
import numpy as np
import os
from astropy.io import fits
import astropy.time as time
import re
import argparse
import pdb

from pydl.pydlutils.yanny import yanny
import yaml
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.gaia import Gaia


from esutil import htm

def old_gproc(mjd=None, observatory=None, outfile=None):
    """Process an MJD of processed guide camera files from plate era into a single file 

    Parameters:
    -----------

    mjd : integer
      MJD (days) of observation
    observatory : str
      which observatory ('apo' or 'lco')

    Comments:
    --------

    Writes gcam-[mjd].fits in local directory with binary table containing
    a row for each guider image.

    Note that 'date-obs' from each guider image header is translated into
    a floating point 'mjd'.
    """

    if observatory == None :
        gdir = os.path.join('/data/gcam',str(mjd))
    elif observatory == 'apo':
        gdir = os.path.join(os.getenv('GCAM_DATA_N'), str(mjd))
    else:
        gdir = os.path.join(os.getenv('GCAM_DATA_S'), str(mjd))

    files = os.listdir(gdir)
    gcam0_dtype = [('gdrms', np.dtype(np.float32)),
                   ('seeing', np.dtype(np.float32)),
                   ('fwhm_median', np.dtype(np.float32)),
                   ('fwhm_mean', np.dtype(np.float32)),
                   ('indx', np.dtype(np.int32)),
                   ('date-obs', np.dtype(np.str_), 21),
                   ('mjd', np.dtype(np.float64)),
                   ('ra', np.dtype(np.float64), 17),
                   ('dec', np.dtype(np.float64), 17),
                   ('xFocal', np.dtype(np.float32), 17),
                   ('yFocal', np.dtype(np.float32), 17),
                   ('focusOffset', np.dtype(np.float32), 17),
                   ('xstar', np.dtype(np.float32), 17),
                   ('ystar', np.dtype(np.float32), 17),
                   ('xCenter', np.dtype(np.float32), 17),
                   ('yCenter', np.dtype(np.float32), 17),
                   ('dx', np.dtype(np.float32), 17),
                   ('dy', np.dtype(np.float32), 17),
                   ('dRA', np.dtype(np.float32), 17),
                   ('dDec', np.dtype(np.float32), 17),
                   ('fwhm', np.dtype(np.float32), 17),
                   ('flux', np.dtype(np.float32), 17),
                   ('mag', np.dtype(np.float32), 17),
                   ('sky', np.dtype(np.float32), 17),
                   ('skymag', np.dtype(np.float32), 17)]
                   
    count = 0
    for file in sorted(files):
        if (file.startswith('proc-')) :
            count = count + 1
    gcam = np.zeros(count, dtype=gcam0_dtype)

    count = 0
    oldname = ''
    ra0=-1
    dec0=-100
    for file in sorted(files):
        #print(file,len(files))
        if (file.startswith('proc-')):
            indx_search = re.search("proc-gimg-([0-9]{4}).fits", file)
            indx = np.int32(indx_search.group(1))
            #header = fitsio.read_header(os.path.join(gdir, file), ext=0)
            hdulist = fits.open(os.path.join(gdir, file))
            header = hdulist[0].header

            if('object' in header['IMAGETYP']):
                gcam['indx'][count] = indx
                try:
                    rescale = (3600. / np.float32(header['PLATSCAL']))
                except:
                    rescale = 0.
                data = hdulist[6].data
                ii = np.nonzero((data['enabled'] == True) &
                                (data['gprobebits'] == 0) &
                                (data['dx'] == data['dx']))[0]
                if len(ii) > 0 :
                    soff = (data['dx'][ii]**2 + data['dy'][ii]**2)
                    rms = np.sqrt(soff.mean())
                    gcam['gdrms'][count] = rms * rescale
                ii = np.nonzero((data['enabled'] == True) &
                                (data['gprobebits'] == 0) &
                                (data['focusOffset'] == 0.) &
                                (data['fwhm'] == data['fwhm']))[0]
                if len(ii) > 0 :
                    fwhm_median = np.median(data['fwhm'][ii])
                    if (fwhm_median == fwhm_median):
                        gcam['fwhm_median'][count] = fwhm_median
                    fwhm_mean = np.mean(data['fwhm'][ii])
                    if (fwhm_mean == fwhm_mean):
                        gcam['fwhm_mean'][count] = fwhm_mean
                for (name, indx) in zip(gcam.dtype.names,
                                        range(len(gcam.dtype.names))):
                    if (indx >= 7):
                        try:
                            gcam[name][count] = data[name]
                        except ValueError:
                            pass
                        except KeyError:
                            pass
                gcam['seeing'][count] = np.float32(header['SEEING'])
                try: dt = str(header['DATE-OBS'].decode())
                except: dt = header['DATE-OBS']
                gcam['date-obs'][count] = dt
                tt = time.Time(dt)
                tt.format = 'mjd'
                gcam['mjd'][count] = tt.value

                if header['NAME'] != oldname :
                    plug,plugheader,sky,stan = getconfig(plugid=header['NAME'],specid=0)
                    guide = np.where(plug['holeType'] == b'GUIDE')[0]
                    oldname = header['NAME']
                gcam['ra'][count,plug['fiberId'][guide]-1] = plug['ra'][guide]
                gcam['dec'][count,plug['fiberId'][guide]-1] = plug['dec'][guide]
                print(gcam['ra'][count],gcam['dec'][count])
                print(gcam['flux'][count])
                if ( np.abs(gcam['ra'][count,0]-ra0)*np.cos(gcam['dec'][count,0]*np.pi/180) > 0.01 or 
                     np.abs(gcam['dec'][count,0]-dec0) > 0.01 ) :
                    ra0=gcam['ra'][count,0]
                    dec0=gcam['dec'][count,0]
                    zp=np.nan
                    z=[]
                    for ra,dec,flux in zip(gcam['ra'][count],gcam['dec'][count],gcam['flux'][count]) :
                      try :
                        z.append(zero(ra,dec,[ra],[dec],[flux],0.002))
                      except: pass
                    print(z)
                    print('new field: ', gcam['ra'][count,0],gcam['dec'][count,0],np.nanmedian(z))
                    pdb.set_trace()

                count = count + 1

    gcam = gcam[0:count - 1]
    if outfile is None: 
        outfile='gcam-{mjd}.fits'.format(mjd=mjd)
    fits.writeto(outfile, gcam, overwrite=True)

    return gcam


def gproc(mjd=None, observatory=None, outfile=None, plot=False):
    """Process an MJD of processed guide camera files from GFA era into a single file

    Parameters:
    -----------

    mjd : integer
      MJD (days) of observation
    instrument : str
      which instrument ('apogee-s' or 'apogee-n')

    Comments:
    --------

    Writes gcam-[mjd].fits in local directory with binary table containing
    a row for each guider image.

    Note that 'date-obs' from each guider image header is translated into
    a floating point 'mjd'.
    """

    if mjd == None :
        mjd = int(time.Time.now().mjd)

    if observatory == None :
        gdir = os.path.join('/data/gcam',str(mjd))
    elif observatory == 'apo':
        gdir = os.path.join(os.getenv('GCAM_DATA_N'), str(mjd))
    else:
        gdir = os.path.join(os.getenv('GCAM_DATA_S'), str(mjd))


    files = glob.glob(gdir+'/proc*')
    
    ids=[]
    for file in files :
        indx_search = re.search(".proc-gimg-gfa.n-([0-9]{4}).fits", file)
        ids.append(indx_search.group(1))

    indexes = sorted(set(ids))

    gcam0_dtype = [
                   ('indx', np.dtype(np.int32)),
                   ('date-obs', np.dtype(np.str_), 21),
                   ('mjd', np.dtype(np.float64)),
                   ('ra', np.dtype(np.float64)),
                   ('dec', np.dtype(np.float64)),
                   ('gdrms', np.dtype(np.float32)),
                   ('dRA', np.dtype(np.float32)),
                   ('dDec', np.dtype(np.float32)),
                   ('dRot', np.dtype(np.float32)),
                   ('dScale', np.dtype(np.float32)),
                   ('corrRA', np.dtype(np.float32)),
                   ('corrDec', np.dtype(np.float32)),
                   ('corrRot', np.dtype(np.float32)),
                   ('corrScale', np.dtype(np.float32)),
                   ('zeropoint', np.dtype(np.float32)),
                   ('nstars', np.dtype(np.int32), 6),
                   ('fwhm', np.dtype(np.float32), 6)
                   ]
    gcam = np.zeros(len(indexes), dtype=gcam0_dtype)
    print(mjd,len(indexes))
    count = 0
    ra0=-1
    dec0=-100
    zp=np.nan
    for index in indexes :
        indx = int(index)
        for cam in range(6) :
          if os.path.exists('{:s}/proc-gimg-gfa{:d}n-{:04d}.fits'.format(gdir,cam+1,int(index))) :
            hdulist = fits.open('{:s}/proc-gimg-gfa{:d}n-{:04d}.fits'.format(gdir,cam+1,int(index)))
            header = hdulist[1].header
            data = hdulist[2].data
            gcam['indx'][count] = int(index)
            try: dt = str(header['DATE-OBS'].decode())
            except: dt = header['DATE-OBS']
            tt = time.Time(dt)
            try :
                gcam['date-obs'][count] = dt
                gcam['mjd'][count] = tt.mjd
                gcam['gdrms'][count] = header['rms']
                gcam['ra'][count] = header['RA']
                gcam['dec'][count] = header['DEC']
                gcam['dRA'][count] = header['DELTARA']
                gcam['dDec'][count] = header['DELTADEC']
                gcam['dRot'][count] = header['DELTAROT']
                gcam['dScale'][count] = header['DELTASCL']
                gcam['corrRA'][count] = header['CORR_RA']
                gcam['corrDec'][count] = header['CORR_DEC']
                gcam['corrRot'][count] = header['CORR_ROT']
                gcam['corrScale'][count] = header['CORR_SCL']
                gcam['fwhm'][count][cam] = header['fwhm']
                gcam['nstars'][count][cam] = len(data)
            except: pass
            if cam == 1 :
                if ( np.abs(gcam['ra'][count]-ra0)*np.cos(gcam['dec'][count]*np.pi/180) > 0.01 or 
                     np.abs(gcam['dec'][count]-dec0) > 0.01 ) :
                  zp=np.nan
                  try:
                    ra0=gcam['ra'][count]
                    dec0=gcam['dec'][count]
                    x=data['x']
                    y=data['y']
                    flux=data['flux']
                    if mjd < 59626 :
                        wcs=WCS(fits.open('{:s}/{:d}/astrometry/proc-gimg-gfa{:d}n-{:04d}.fits.wcs'.format(os.environ['GCAM_DATA_N'],mjd,cam+1,indx))[0].header)
                    else :
                        wcs=WCS(fits.open('{:s}/{:d}/astrometry/gimg-gfa{:d}n-{:04d}.wcs'.format(os.environ['GCAM_DATA_N'],mjd,cam+1,indx))[0].header)
                    out=wcs.pixel_to_world(x,y)
                    ra= out[:].ra.value
                    dec=out[:].dec.value 
                    zp=zero(header['CRVAL1'],header['CRVAL2'],ra,dec,flux,0.1)
                    print('new field: ', gcam['ra'][count],gcam['dec'][count],zp)
                  except: pass
            gcam['zeropoint'][count] = zp
        count += 1

    gcam = gcam[0:count - 1]
    if outfile is None: 
        outfile=os.path.join(gdir,'gcam-{mjd}.fits'.format(mjd=mjd))
    fits.writeto(outfile, gcam, overwrite=True)

    if plot :
        fig,ax=plots.multi(3,2,figsize=(12,8))
        plots.plotc(ax[0,0],gcam['mjd'],gcam['gdrms'],np.median(gcam['fwhm'],axis=1),yr=[0,2],xt='MJD',yt='gdrms')
        for i in range(6) : plots.plotp(ax[0,1],gcam['mjd'],gcam['fwhm'][:,i],label='Camera {:d}'.format(i+1),color=None,xt='MJD',yt='fwhm',yr=[0,5])
        ax[0,1].legend(fontsize='xx-small')
        plots.plotc(ax[1,0],gcam['mjd'],gcam['dRA'],np.median(gcam['fwhm'],axis=1),yr=[0,2],xt='MJD',yt='dRA')
        plots.plotc(ax[1,1],gcam['mjd'],gcam['dDec'],np.median(gcam['fwhm'],axis=1),yr=[0,2],xt='MJD',yt='dDec')
        plots.plotc(ax[1,2],gcam['mjd'],gcam['dScale'],np.median(gcam['fwhm'],axis=1),yr=[0.9995,1.0005],xt='MJD',yt='dScale')
        ax[1,2].get_yaxis().get_major_formatter().set_useOffset(False)
        fig.tight_layout()
        fig.savefig(outfile.replace('.fits','.png'))

    return gcam

def zero(ra0,dec0,ra,dec,flux,rad) :

    gaia_posn = get_gaia(ra0,dec0,rad)
    h=htm.HTM()
    maxrad=3./3600.
    gmag=gaia_posn['phot_g_mean_mag']
    inst=-2.5*np.log10(flux)
    m1,m2,rad=h.match(ra,dec, gaia_posn['ra'], gaia_posn['dec'], maxrad,maxmatch=1) 
    gd = np.where((gmag[m2]>10)&(gmag[m2]<17))[0]
    if len(gd) < 1 : 
        print('no good stars: ', len(m2), len(gd))
        return np.nan
    print(gmag)
    print(flux)

    return np.nanmedian(gmag[m2[gd]]-inst[m1[gd]])

def guider_zero(gcams,mjd,mjd_obs,exptime) :

  # find appropriate guider exposures for this exposure
  gnearest = np.argmin(np.abs(gcams['mjd']-mjd_obs))
  gexp = np.where((gcams['mjd']>mjd_obs) &
                  (gcams['mjd']<mjd_obs+exptime/3600./24.) )[0]
  if len(gexp) == 0 : gexp = np.array([gnearest])

  # loop through the exposures
  gzero = []
  secz = []
  #for indx in indxs :
  for gcam in gcams[gexp] :

    if mjd < 59500 :
        # do GAIA query for first image only
        gfa=0
        outfile='through/{:d}/proc-gimg-{:04d}.fits'.format(mjd,gcams[0]['indx'])
        try: os.mkdir('through/{:d}'.format(mjd))
        except FileExistsError : pass
        ra = gcam['ra']
        dec = gcam['dec']
        flux = gcam['flux']
        try :
            gaia_posn = Table.read(outfile)
        except :
            out=[]
            for r,d in zip(ra,dec) :
                out.append(get_gaia(r,d,15/3600.))
            gaia_posn = vstack(out)
            gaia_posn['ra','dec','phot_g_mean_mag'].write(outfile)

    else :
        gfa = 2
        indx = gcam['indx']
        infile='{:s}/{:d}/proc-gimg-gfa{:d}n-{:04d}.fits'.format(os.environ['GCAM_DATA_N'],mjd,gfa,indx)
        try : a=fits.open(infile)
        except : continue
 
        try : 
            if mjd < 59626 :
                wcs=WCS(fits.open('{:s}/{:d}/astrometry/proc-gimg-gfa{:d}n-{:04d}.fits.wcs'.format(os.environ['GCAM_DATA_N'],mjd,gfa,indx))[0].header)
            else :
                wcs=WCS(fits.open('{:s}/{:d}/astrometry/gimg-gfa{:d}n-{:04d}.wcs'.format(os.environ['GCAM_DATA_N'],mjd,gfa,indx))[0].header)
            out=wcs.pixel_to_world(a[2].data['x'],a[2].data['y'])
            ra= out[:].ra.value
            dec=out[:].dec.value 
        except :
            print('error with WCS transformation...{:s}/{:d}/astrometry/proc-gimg-gfa{:d}n-{:04d}.fits.wcs'.format(os.environ['GCAM_DATA_N'],mjd,gfa,indx))
            continue
        ra0 = a[1].header['CRVAL1']
        dec0 = a[1].header['CRVAL2']
        rad = 0.1
        flux=a[2].data['flux']

        # do GAIA query for first image only
        outfile='through/{:d}/proc-gimg-gfa{:d}n-{:04d}.fits'.format(mjd,gfa,gcams[0]['indx'])
        try: os.mkdir('through/{:d}'.format(mjd))
        except FileExistsError : pass
        try :
            gaia_posn = Table.read(outfile)
        except :
            gaia_posn = get_gaia(ra0,dec0,rad)
            gaia_posn['ra','dec','phot_g_mean_mag'].write(outfile)
        
    h=htm.HTM()
    maxrad=3./3600.
    m1,m2,rad=h.match(ra,dec, gaia_posn['ra'], gaia_posn['dec'], maxrad,maxmatch=1) 
    #plt.figure()
    #plt.plot(ra,dec,'ro')
    #plt.plot(gaia_posn['ra'],gaia_posn['dec'],'go')
    #plt.draw()
    #print(a[1].header['CRPIX1'],a[1].header['CRPIX2'])
    #pdb.set_trace()
    #plt.close()
    gmag=gaia_posn['phot_g_mean_mag']
    inst=-2.5*np.log10(flux)
    gd = np.where((gmag[m2]>13)&(gmag[m2]<17))[0]
    #print(infile,len(m2),len(gd))
    if len(gd) < 1 : 
        print('no good stars: ', len(m2), len(gd))
        continue

    gzero.append(np.nanmedian(gmag[m2[gd]]-inst[m1[gd]]))
    #secz.append(1./np.cos((90-a[1].header['ALT'])*np.pi/180.))

  print('{:d} {:8.2f} {:8.2f}'.format(gfa,np.median(gzero),np.median(secz)))

  if len(gzero) == 0 : return (20.,30.,0.),0.
  elif len(gzero) == 1 : return((gzero[0],gzero[0]-5,gzero[0]+5),np.median(secz))
  else : return np.percentile(gzero,[50,25,75])    #,np.median(secz)


def get_gaia(ra,dec,rad) :

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    radius = u.Quantity(rad, u.deg)
    print('radius: ', radius)
    Gaia.ROW_LIMIT = -1
    job = Gaia.cone_search_async(coord, radius)
    r = job.get_results()
    return r



# Process guide camera images into a single roll-up file with statistics.

def getconfig(config_id=None,plugid=None,specid=1) :
    """ read confSummary or plPlugMap file, return data
    """
    if config_id is not None :
        plug,header=config(config_id,specid=specid,useconfF=True,useparent=False)
        sky=np.where(plug['category'] == b'sky_boss')[0]
        stan=np.where(np.char.find(plug['category'].astype(str),'standard') >= 0)[0]
        if specid == 1 :
            # substitude GAIA transformed mags for gri for gaia_g < 15
            x = plug['bp_mag']-plug['rp_mag']
            x2 = x * x
            x3 = x * x * x
            gaia_G = plug['gaia_g_mag']
            j = np.where(gaia_G < 15)[0]
            plug['mag'][j,1] = -1 * (0.13518 - 0.46245 * x[j] - 0.25171 * x2[j] + 0.021349 * x3[j]) + gaia_G[j]
            plug['mag'][j,2] = -1 * (-0.12879 + 0.24662 * x[j] - 0.027464 * x2[j] - 0.049465 * x3[j]) + gaia_G[j]
            plug['mag'][j,3] = -1 * (-0.29676 + 0.64728 * x[j] - 0.10141 * x2[j]) + gaia_G[j]
    elif plugid is not None :
        plug,header=config(os.environ['MAPPER_DATA_N']+'/'+plugid.split('-')[1]+'/plPlugMapM-'+plugid+'.par',specid=specid,struct='PLUGMAPOBJ')
        sky=np.where(plug['objType'] == b'SKY')[0]
        stan=np.where(np.char.find(plug['objType'].astype(str),'STD') >= 0)[0]
        plug=Table(plug)
        plug['h_mag']=np.nan

        # get plateHoles file
        plate=int(plugid.split('-')[0])
        # substitute H mag from plateHoles
        holes=yanny('{:s}/plates/{:04d}XX/{:06d}/plateHolesSorted-{:06d}.par'.format(
                  os.environ['PLATELIST_DIR'],plate//100,plate,plate))
        h=htm.HTM()
        m1,m2,rad=h.match(plug['ra'],plug['dec'],
                  holes['STRUCT1']['target_ra'],holes['STRUCT1']['target_dec'],
                  0.1/3600.,maxmatch=500)
        if specid == 1 :
            if int(header['plateId']) >= 15000 :
                corr = fits.open(
                         os.environ['IDLSPEC2D_DIR']+'/catfiles/Corrected_values_plate{:s}_design{:s}.fits'.format(
                         header['plateId'],header['designid']))[1].data
                h1,h2,rad=h.match(plug['ra'][m1],plug['dec'][m1], corr['RA'],corr['DEC'],
                          0.1/3600.,maxmatch=500)
                bd=np.where(holes['STRUCT1']['catalogid'][m2[h1]] == np.array(corr['Original_CatalogID'],dtype=np.int64)[h2])[0]
                j=np.where(corr[h2[bd]]['Mag_Change'])[0]
                print(plugid,len(bd),len(j))
                plug['mag'][m1[h1[bd]],1] = corr['gmag'][h2[bd]]
                plug['mag'][m1[h1[bd]],2] = corr['rmag'][h2[bd]]
                plug['mag'][m1[h1[bd]],3] = corr['imag'][h2[bd]]
                plug['mag'][m1[h1[bd]],4] = corr['zmag'][h2[bd]]
        elif specid == 2 :
            plug['h_mag'][m1] = holes['STRUCT1']['tmass_h'][m2]
    else :
        raise ValueError('either config_id or plugid needs to be set')

    return np.array(plug),header,sky,stan

def config(cid,specid=2,struct='FIBERMAP',useparent=True,useconfF=False) :
    """ Get FIBERMAP structure from configuration file for specified instrument
           including getting parent_configuration if needed (for scrambled configurations)
    """
    if useconfF : confname='confSummaryF'
    else :confname='confSummary'
    if isinstance(cid,str):
        conf = yanny(cid)
    else :
        conf = yanny(os.environ['SDSSCORE_DIR']+'/apo/summary_files/{:04d}XX/{:s}-{:d}.par'.format(cid//100,confname,cid))
        if useparent :
            try :
                parent = int(conf['parent_configuration'])
                if parent > 0 :
                    conf = yanny(os.environ['SDSSCORE_DIR']+'/apo/summary_files/{:04d}XX/{:s}-{:d}.par'.format(parent//100,confname,parent))
            except :  pass

    if conf == None or len(conf) == 0 :
        raise FileNotFoundError('error opening file',cid)

    gd =np.where((conf[struct]['spectrographId'] == specid) & (conf[struct]['fiberId'] > 0) )[0]
    return conf[struct][gd],conf.new_dict_from_pairs()


