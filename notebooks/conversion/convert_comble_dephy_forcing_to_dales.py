#!/bin/python3
# =============================================
#  COMBLE case conversion: COMBLE -> Dales tb
#  script for the conversion of input files
# =============================================

# ----- basic setup ----------------
# input file
dephy_filename = '../forcing/COMBLE_INTERCOMPARISON_FORCING_V2.3.nc'

# output directory
out_dir = '../../dales/'
parent_dir = './'

# output scm file 
out_scmname = 'scm_comble_v2.3_dales_v002.nc'
out_profname=     'prof.inp.001'

# case info
caseatt = {'date':"20200312",
           'source':"comble-mip",
           'source_datafiles_prof':'https://github.com/ARM-Development/comble-mip/tree/main/data_files',
           'source_datafiles_sfc':'https://github.com/ARM-Development/comble-mip/tree/main/data_files',
           'creator':"Jan Chylik and Roel Neggers, IGMK, University of Cologne"
            }


# whether to write outpus 
is_verbose = True


# ------ library imports ---------------

# import libraries 
# later turn into module 
import numpy              as np
import matplotlib.pyplot  as plt
import os 
import datetime           as dt
from   netCDF4        import Dataset

# append path to import the module

# ----- inner setup ----------------

comsetup = {}
comsetup ['repeatdim'] = (0,'t0','time') 

def_fill_value = -999.

prof4_com = 'zw_grid' 


vari_dales={}
# vari_dales['init']=['height','thetal','qt','u','v','tke_init']
vari_dales['init']=['height','thetal','qt','u','v','tke']
vari_dales['defval']=[0,273.0,0.5e-5,1.0,1.0, 0.1]

# formatting thresholds
val_fill = 0.0
max_const = 1.0e5
med_const = 1.6e4
min_const = 1.0e-3
eps_const = 1.0e-13
#
len_pos = 14        # width of the text column
len_dec_std = 5      # standart number of dec places
len_dec_exp = 4      # number of decimal places for exp format 



# ------ future module --------------


# function definitions
def addnlevp1(varname,nlv):
    return nlv+1

def nlevels(varname):
    return np.array(range(1,1+tbnc[varname].shape[0]),dtype=np.int32)

def toint (varname,val):
    return np.int32(val)

def convtime (varname,val):
    # add offset form 00 of the day
    time1day =val #dt1day+val
    return np.int32(time1day)

def convepoch (varname,val):
    # add epoch offset 
    time1970day = dtstamp+val
    return np.int32(time1970day)


def convdateint (varname,val):
    # as a list, because array conversion does not work
    tlist = []
    # I know it is terrible, but utcfromtimestamp does not like arrays
    for j in range(val.shape[0]):
        tlist.append(int(dt.datetime.utcfromtimestamp(dtstamp+val[j]).strftime("%Y%m%d")))
    # add epoch offset 
    return np.array(tlist,dtype=np.int32)


def repeatfill(varname,val):
    arr = val*np.ones(tbnc[varname].shape)
    return arr    


def filldef(varname):
    arr = tbvars[varname]['def_value']*np.ones(tbnc[varname].shape)
    return arr



def write_inpfile (fileout, val, uptext='! \n ! \n', mitext='! \n ! \n'):    
    # writes inpfile 
    n_vari = val.shape[1]  # length of the list of variables
    n_levs = val.shape[0]
    #
    #later adjust that it works for time dependent variables as well
    #
    # -------------------------
    # prepare array to write 
    # -------------------------
    # 
    # -------------------------
    # prepare format 
    # -------------------------
    # 2. preparing a formatting string
    form_str = ""
    for jvar in range(0,n_vari):
        # - look up the max and min of the variable
        # minval = np.amin(np.abs(val_out[:,jvar]))
        maxval = np.amax(np.abs(val[:,jvar]))
        #
        if( (maxval>max_const) ):
            # - exponential format for high numbers
            str_part = '%s%d.%de ' % ('%',len_pos,len_dec_exp)
        elif ( (maxval<min_const) & (maxval> eps_const)):
            # - exponential format for low numbers
            str_part = '%s%d.%de ' % ('%',len_pos,len_dec_exp) 
        elif ( (maxval>min_const) & (maxval<med_const) ):
            str_part = '%s%d.%df ' % ('%',len_pos,len_dec_exp)
        else:        
            # - standard format for lower numbers
            str_part = '%s%d.%df ' % ('%',len_pos,len_dec_std) 
        #
        # put it together
        form_str += str_part
    # -------------------------
    # writing the file 
    # -------------------------           
    # prepare first two lines:
    # open file to write into 
    ofile = open(fileout, 'wb')
    # write it 
    ofile.write(uptext.encode()) # so it can be writtne in block mode
    # writing the arrays with numpy.savetxt
    np.savetxt( ofile, val, form_str )
    #
    # close the file
    ofile.close()
    statsig = 0
    return statsig


def write_prof_inp (zvec,fileout='prof.inp.001'):
    # prepare dales prof.inp.iexpr file
    # 
    # number of levels
    n_levs = np.shape(zvec)[0]
    # start writing output
    twolines_prof = "# Input profiles, %s, Nz = %d  \n" %('title' , n_levs)  
    twolines_prof += "#  "+"    ".join(vari_dales['init'])+" \n"
    val_out = np.zeros((n_levs,len(vari_dales['init'])), dtype = 'float')
    # first column are altitudes
    jvar = 0
    val_out[:,jvar] = zvec[:]
    for jvar in range(1,val_out.shape[1]):
        val_out[:,jvar] = vari_dales['defval'][jvar]*np.ones((n_levs))
    write_inpfile (fileout,val_out ,uptext=twolines_prof)
    statsig = 0
    # 
    return statsig 


tb_dims = {}
tb_dims['time']="time"
tb_dims['nlev']="vertical levels"
tb_dims['nlev1']="vertical half-levels"
tb_dims['nlevs']="vertical levels in soil model"
tb_dims['nsv']="scalar variables"
tb_dims['nDS']="number of dropsondes"


tb_fun = {}
# dictionary of functions
tb_fun['toint']      = toint
tb_fun['convtime']   = convtime
tb_fun['convepoch']  = convepoch
tb_fun['convdateint']= convdateint
tb_fun['addnlevp1']  = addnlevp1
tb_fun['nlevels']    = nlevels
tb_fun['repeatfill'] = repeatfill
tb_fun['filldef']    =  filldef



tb4com_dims = {}
#  [dales variable] = (conversion type,key or function, input keys)
tb4com_dims['time']    = ('dim','time')
tb4com_dims['nlev']    = ('dim','lev')
tb4com_dims['nlevp1']  = ('fun','addnlevp1',('lev',))  # double check the function call works
# tb4com_dims['nlevs']   = ('dim','time')



tb4com_vars = {}
#  [dales variable] = (conversion type,key or function, input keys)
tb4com_vars['time']        = ('fun','convtime',('time',))
tb4com_vars['base_time']   = ('fun','convepoch',('time',))
tb4com_vars['date']        = ('fun','convdateint',('time',))
tb4com_vars['second']      = ('fun','toint',('time',))
tb4com_vars['nlev']        = ('fun','nlevels',())
tb4com_vars['nlevp1']      = ('fun','nlevels',())
#tb4com_vars['nlevs']       = ('fun','nlevels',())
tb4com_vars['lat']         = ('fun','repeatfill',('lat',))
tb4com_vars['lon']         = ('fun','repeatfill',('lon',))
tb4com_vars['height_f']    = ('var','lev')
tb4com_vars['pressure_f']  = ('var','pressure')
# reference values
tb4com_vars['lat_grid']         = ('fun','repeatfill',('lat_ref',))
tb4com_vars['lon_grid']         = ('fun','repeatfill',('lon_ref',))
# thermodynamic variables
tb4com_vars['u']           = ('var','u')
tb4com_vars['v']           = ('var','v')
tb4com_vars['t']           = ('var','temp')
tb4com_vars['q']           = ('var','qv')
tb4com_vars['ql']          = ('fun','filldef',())  # ('var','ql')
tb4com_vars['qi']          = ('fun','filldef',())  # ('var','qi')
tb4com_vars['cloud_fraction']= ('fun','filldef',())
tb4com_vars['omega']        = ('fun','filldef',())  # ('var','wap')
tb4com_vars['o3']           = ('var','o3')
tb4com_vars['t_local']       = ('var','temp')
tb4com_vars['q_local']       = ('var','qv')
tb4com_vars['ql_local']      = ('fun','filldef',()) # ('var','ql')
tb4com_vars['qi_local']      = ('fun','filldef',()) # ('var','qi')
tb4com_vars['u_local']       = ('var','u')
tb4com_vars['v_local']       = ('var','v')
tb4com_vars['cc_local']      = ('fun','filldef',())
# tendencies
tb4com_vars['tadv']           = ('fun','filldef',())  # ('var','tnta_adv')
tb4com_vars['qadv']           = ('fun','filldef',())  # ('var','tnqa_adv')
tb4com_vars['uadv']           = ('fun','filldef',())  # ('var','tnua_adv')
tb4com_vars['vadv']           = ('fun','filldef',())  # ('var','tnva_adv')
tb4com_vars['ladv']           = ('fun','filldef',())
tb4com_vars['iadv']           = ('fun','filldef',())
tb4com_vars['aadv']           = ('fun','filldef',())
# geostrophic wind
tb4com_vars['ug']           = ('var','ug')
tb4com_vars['vg']           = ('var','vg')
# surface variables
tb4com_vars['ps']           = ('var','ps')
tb4com_vars['fradSWnet']    = ('fun','filldef',())
tb4com_vars['fradLWnet']    = ('fun','filldef',())
tb4com_vars['albedo']       = ('fun','filldef',())
tb4com_vars['sfc_sens_flx']  = ('fun','filldef',())
tb4com_vars['sfc_lat_flx']   = ('fun','filldef',())
tb4com_vars['mom_rough']    = ('att','z0')
tb4com_vars['heat_rough']   = ('att','z0h')
tb4com_vars['t_skin']       = ('var','ts')

## optional variables
#x lat_grid
#x lon_grid
# pressure_h
# gz_f
# 


# testbed variables
# output generated by: 
#tbfile = '/work/jchylik/cases/mosaic_aca/mosaic-aca/mosaic_aca_20200911_def/scm_in.ERA5_CDS_MOSAiC_coor-LagrangianLagranto-traj0-p950_domain2.0x2.0_20200910_ndays4_cut-014hr_surf_qcor4adv0.nc'
#tbinnc = Dataset(tbfile,'r')
#for jvar in tbinnc.variables.keys():
    #vartype = tbinnc[jvar].dtype          # get format
    #dims = tbinnc[jvar].dimensions        # get dimensions
    #attlist=[]
    #for attr in tbinnc[jvar].ncattrs():   # loop over attributes
        #attlist.append("'{}':'{}'".format(attr,tbinnc[jvar].getncattr(attr)))
    #print("tbvars['{}'] = {{ 'dtype':'{}', 'dimensions': {},{} }}".format(jvar,vartype,dims,",".join(attlist)))
#tbinnc.close()

atts_tbvars = ['units','long_name']

tbvars ={}
tbvars['time'] = { 'dtype':'int32', 'dimensions': ('time',),'long_name':'time','units':'seconds since 00 UTC on first day','_FillValue':-999 }
tbvars['base_time'] = { 'dtype':'int32', 'dimensions': ('time',),'long_name':'epoch time','units':'seconds since 1-1-1970 00:00','_FillValue':-999 }
tbvars['date'] = { 'dtype':'int32', 'dimensions': ('time',),'long_name':'date','units':'yyyymmdd','_FillValue':-999 }
tbvars['second'] = { 'dtype':'int32', 'dimensions': ('time',),'long_name':'seconds since start of sequence','units':'s','_FillValue':-999 }
tbvars['lat'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'latitude','units':'degrees North','_FillValue':-999.0 }
tbvars['lon'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'longitude','units':'degrees East','_FillValue':-999.0 }
tbvars['lat_grid'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'latitude of closest IFS gridpoint','units':'degrees North','_FillValue':-999.0 }
tbvars['lon_grid'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'longitude of closest IFS gridpoint','units':'degrees East','_FillValue':-999.0 }
tbvars['height_f'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'full level height','units':'m','_FillValue':-999.0 }
tbvars['nlev'] = { 'dtype':'int32', 'dimensions': ('nlev',),'long_name':'model full levels','units':'positive integers','_FillValue':-999 }
tbvars['height_h'] = { 'dtype':'float32', 'dimensions': ('time', 'nlevp1'),'long_name':'half level height','units':'m','_FillValue':-999.0 }
tbvars['nlevp1'] = { 'dtype':'int32', 'dimensions': ('nlevp1',),'long_name':'model half levels','units':'positive integers','_FillValue':-999 }
tbvars['ps'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'surface pressure','units':'Pa','_FillValue':-999.0 }
tbvars['pressure_f'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'full level pressure','units':'Pa','_FillValue':-999.0 }
tbvars['pressure_h'] = { 'dtype':'float32', 'dimensions': ('time', 'nlevp1'),'long_name':'half level pressure','units':'Pa','_FillValue':-999.0 }
tbvars['gz_f'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'geoptential height','units':'m2/s2','_FillValue':-999.0 }
tbvars['u'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'zonal wind (domain averaged)','units':'m/s','_FillValue':-999.0 }
tbvars['v'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'meridional wind (domain averaged)','units':'m/s','_FillValue':-999.0 }
tbvars['t'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'temperature (domain averaged)','units':'K','_FillValue':-999.0 }
tbvars['q'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'water vapor mixing ratio (domain averaged)','units':'kg/kg','_FillValue':-999.0 }
tbvars['ql'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'liquid water mixing ratio (domain averaged)','units':'kg/kg', 'def_value':0.0,'_FillValue':-999.0 }
tbvars['qi'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'ice water mixing ratio (domain averaged)','units':'kg/kg', 'def_value':0.0,'_FillValue':-999.0 }
tbvars['cloud_fraction'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'cloud fraction (domain averaged)','units':'0-1','_FillValue':-999.0, 'def_value':0.0}
tbvars['omega'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'large-scale pressure velocity (domain averaged)','units':'Pa/s','def_value':0.0,'_FillValue':-999.0 }
tbvars['o3'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'ozone mass mixing ratio (domain averaged)','units':'kg/kg','_FillValue':-999.0 }
tbvars['t_local'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'temperature (at domain midpoint)','units':'K','_FillValue':-999.0 }
tbvars['q_local'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'water vapor specific humidity (at domain midpoint)','units':'kg/kg','_FillValue':-999.0 }
tbvars['ql_local'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'liquid water specific humidity (at domain midpoint)','units':'kg/kg', 'def_value':0.0,'_FillValue':-999.0 }
tbvars['qi_local'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'ice water specific humidity (at domain midpoint)','units':'kg/kg', 'def_value':0.0,'_FillValue':-999.0 }
tbvars['u_local'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'zonal wind (at domain midpoint)','units':'m/s','_FillValue':-999.0 }
tbvars['v_local'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'meridional wind (at domain midpoint)','units':'m/s','_FillValue':-999.0 }
tbvars['cc_local'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'long_name':'cloud fraction (at domain midpoint)','units':'0-1','_FillValue':-999.0, 'def_value':0.0 }
tbvars['tadv'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'lagrangian':'Lagrangian setup: horizontal advection calculated using velocity relative to wind on trajectory (u_traj,v_traj)','info':'derived at pressure levels','long_name':'tendency in temperature due to large-scale horizontal advection','units':'K/s','_FillValue':-999.0, 'def_value':0.0 }
tbvars['qadv'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'lagrangian':'Lagrangian setup: horizontal advection calculated using velocity relative to wind on trajectory (u_traj,v_traj)','info':'derived at pressure levels','long_name':'tendency in water vapor due to large-scale horizontal advection','units':'kg/kg/s','_FillValue':-999.0, 'def_value':0.0 }
tbvars['uadv'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'lagrangian':'Lagrangian setup: horizontal advection calculated using velocity relative to wind on trajectory (u_traj,v_traj)','info':'derived at pressure levels','long_name':'tendency in zonal wind due to large-scale horizontal advection','units':'m/s2','_FillValue':-999.0 , 'def_value':0.0}
tbvars['vadv'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'lagrangian':'Lagrangian setup: horizontal advection calculated using velocity relative to wind on trajectory (u_traj,v_traj)','info':'derived at pressure levels','long_name':'tendency in meridional wind due to large-scale horizontal advection','units':'m/s2','_FillValue':-999.0, 'def_value':0.0 }
tbvars['ug'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'interpolation':'above 5 hPa the geostrophic wind is equal to the real wind','info':'derived at pressure levels','long_name':'geostrophic wind - zonal component','units':'m/s','_FillValue':-999.0, 'def_value':0.0 }
tbvars['vg'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'interpolation':'above 5 hPa the geostrophic wind is equal to the real wind','info':'derived at pressure levels','long_name':'geostrophic wind -meridional component','units':'m/s','_FillValue':-999.0, 'def_value':0.0 }
tbvars['ladv'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'lagrangian':'Lagrangian setup: horizontal advection calculated using velocity relative to wind on trajectory (u_traj,v_traj)','long_name':'tendency in liquid water spec hum due to large-scale horizontal advection','units':'kg/kg/s','_FillValue':-999.0, 'def_value':0.0 }
tbvars['iadv'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'lagrangian':'Lagrangian setup: horizontal advection calculated using velocity relative to wind on trajectory (u_traj,v_traj)','long_name':'tendency in frozen water due to large-scale horizontal advection','units':'kg/kg/s','_FillValue':-999.0 , 'def_value':0.0}
tbvars['aadv'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'lagrangian':'Lagrangian setup: horizontal advection calculated using velocity relative to wind on trajectory (u_traj,v_traj)','long_name':'tendency in cloud fraction due to large-scale horizontal advection','units':'1/s','_FillValue':-999.0 , 'def_value':0.0}
tbvars['timDS'] = { 'dtype':'int32', 'dimensions': ('nDS',),'units':'seconds since 1-1-1970 00:00','long_name':'time at trajectory reference point','info':'the reference point is the space-time coordinate from which the trajectory is calculated' }
tbvars['latDS'] = { 'dtype':'float32', 'dimensions': ('nDS',),'units':'degrees North','long_name':'latitude at trajectory reference point','info':'the reference point is the space-time coordinate from which the trajectory is calculated' }
tbvars['lonDS'] = { 'dtype':'float32', 'dimensions': ('nDS',),'units':'degrees East','long_name':'longitude at trajectory reference point','info':'the reference point is the space-time coordinate from which the trajectory is calculated' }
tbvars['time_traj'] = { 'dtype':'int32', 'dimensions': ('time',),'long_name':'time values at trajectory waypoints','units':'seconds since 1-1-1970 00:00','_FillValue':-999 }
tbvars['lat_traj'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'latitude of trajectory waypoints','units':'degrees North','_FillValue':-999.0 }
tbvars['lon_traj'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'longitude of trajectory waypoints','units':'degrees East','_FillValue':-999.0 }
tbvars['p_traj'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'pressure level at which trajectory was calculated','units':'hPa','_FillValue':-999.0 }
tbvars['u_traj'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'zonal wind at trajectory waypoints','units':'m/s','_FillValue':-999.0 }
tbvars['v_traj'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'meridional wind at trajectory waypoints','units':'m/s','_FillValue':-999.0 }
tbvars['fradSWnet'] = { 'dtype':'float32', 'dimensions': ('time', 'nlevp1'),'long_name':'radiative flux - net short wave','units':'W/m2','_FillValue':-999.0, 'def_value':0.0 }
tbvars['fradLWnet'] = { 'dtype':'float32', 'dimensions': ('time', 'nlevp1'),'long_name':'radiative flux - net long wave','units':'W/m2','_FillValue':-999.0, 'def_value':0.0 }
tbvars['albedo'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'albedo','units':'0-1','_FillValue':-999.0, 'def_value':0.0 }
tbvars['mom_rough'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'roughness length for momentum','units':'m','_FillValue':-999.0 }
tbvars['heat_rough'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'roughness length for heat','units':'m','_FillValue':-999.0 }
tbvars['t_skin'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'skin temperature','units':'K','_FillValue':-999.0 }
tbvars['q_skin'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'m of water','units':'skin reservoir content','_FillValue':-999.0 }
tbvars['high_veg_type'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'high vegetation type','units':'-','_FillValue':-999.0 }
tbvars['low_veg_type'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'low vegetation type','units':'-','_FillValue':-999.0 }
tbvars['high_veg_cover'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'high vegetation cover','units':'0-1','_FillValue':-999.0 }
tbvars['low_veg_cover'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'low vegetation cover','units':'0-1','_FillValue':-999.0 }
tbvars['high_veg_lai'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'leaf area index of high vegetation','units':'-','_FillValue':-999.0 }
tbvars['low_veg_lai'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'leaf area index of low vegetation','units':'-','_FillValue':-999.0 }
tbvars['snow'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'snow depth','units':'m, liquid equivalent','_FillValue':-999.0 }
tbvars['t_snow'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'snow temperature','units':'K','_FillValue':-999.0 }
tbvars['albedo_snow'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'snow albedo','units':'0-1','_FillValue':-999.0 }
tbvars['density_snow'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'snow density','units':'kg/m3', 'def_value':0.0,'_FillValue':-999.0 }
tbvars['sfc_sens_flx'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'surface sensible heat flux','units':'W/m2', 'def_value':0.0,'_FillValue':-999.0 }
tbvars['sfc_lat_flx'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'surface latent heat flux','units':'W/m2', 'def_value':0.0,'_FillValue':-999.0 }
tbvars['h_soil'] = { 'dtype':'float32', 'dimensions': ('nlevs',),'_FillValue':-999.0,'units':'m','long_name':'soil layer thickness' }
tbvars['nlevs'] = { 'dtype':'int32', 'dimensions': ('nlevs',),'_FillValue':-999,'long_name':'soil levels' }
tbvars['t_soil'] = { 'dtype':'float32', 'dimensions': ('time', 'nlevs'),'long_name':'soil layer temperature','units':'K','_FillValue':-999.0 }
tbvars['q_soil'] = { 'dtype':'float32', 'dimensions': ('time', 'nlevs'),'long_name':'soil moisture','units':'m3/m3','_FillValue':-999.0 }
tbvars['lsm'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'land sea mask','units':'-','_FillValue':-999.0 }
tbvars['sea_ice_frct'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'sea ice fraction','units':'0-1','_FillValue':-999.0 }
tbvars['t_sea_ice'] = { 'dtype':'float32', 'dimensions': ('time', 'nlevs'),'long_name':'sea ice temperature','units':'K','_FillValue':-999.0 }
tbvars['open_sst'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'open sea surface temperature','units':'K','_FillValue':-999.0 }
tbvars['sdor'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'subgrid-scale orography - standard deviation','units':'m2/s2','_FillValue':-999.0 }
tbvars['isor'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'subgrid-scale orography - anisotropy','units':'01-','_FillValue':-999.0 }
tbvars['anor'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'subgrid-scale orography - orientation/angle of steepest gradient','units':'degree','_FillValue':-999.0 }
tbvars['slor'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'subgrid-scale orography - mean slope','units':'m/m','_FillValue':-999.0 }
tbvars['orog'] = { 'dtype':'float32', 'dimensions': ('time',),'long_name':'orography - surface geopotential','units':'m2/s2','_FillValue':-999.0 }
tbvars['sv'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev', 'nsv'),'long_name':'tracers','units':'whatever','_FillValue':-999.0 }
tbvars['t_skin_seaice'] = { 'dtype':'float32', 'dimensions': ('time',),'units':'K','long_name':'skin temperature - sea ice (scaled to match with MetCity/ASFS observations)' }
tbvars['t_skin_ocean'] = { 'dtype':'float32', 'dimensions': ('time',),'units':'K','long_name':'skin temperature - ocean' }
tbvars['n_ccn'] = { 'dtype':'float32', 'dimensions': ('time', 'nlev'),'units':'/m3','long_name':'CCN concentration (initial)' }




# ----- dephy-tb setup ----------------
#-> move to the mpodule file



# -> later add importing dephydales module


# ------ preparation  ---------------

# check drive file
if os.path.exists(dephy_filename):
    pass
else:
    print("ERROR: Driver file dephy_filename={} not found".format(dephy_filename))
    

# prepare out directory
if os.path.exists(out_dir):
    print("Directory out_dir={} for Dales files already exists".format(out_dir))
else:
    os.mkdir(out_dir)
    if os.path.exists(out_dir):
        print("Creating directory out_dir={} for Dales files".format(out_dir))
    else:
        print("ERROR: Cannot create directory out_dir={} for Dales files".format(out_dir))


# prepare output file 
tbnc = Dataset(out_dir+out_scmname,'w')

# open the input file
dnc = Dataset(dephy_filename,'r') 


# get starting time 
datestringin = dnc.getncattr('startDate')
#nowdt=dt.datetime.now()
#nowdt.strftime('%Y-%m-%d %H:%M:%S')
dtin = dt.datetime.strptime( datestringin, '%Y-%m-%d %H:%M:%S')

# get day of the year
dayoftheyear = dtin.strftime('%j')

# fixing the offset because tzinfo DOES NOT work CORRECTLY
# 1. get computational time offset
ctimeoffset = dt.datetime.fromtimestamp(dtin.timestamp()) - dt.datetime.utcfromtimestamp(dtin.timestamp())
# 2. correct the offset
dtcase = dtin+ctimeoffset
# seconds since 1-1-1970 00:00
dtstamp = dt.datetime.timestamp(dtcase)


# UTC on first day 
dt1day = 3600*dtin.hour+60*dtin.minute






# set expected testbed attributes
for attr in caseatt.keys():
    valattr = caseatt[attr]
    tbnc.setncattr(attr,valattr)
    if is_verbose:
        print("{} = {}".format(attr,valattr))


# copy global attributes
for attr in  dnc.ncattrs():
    valattr= dnc.getncattr(attr)
    tbnc.setncattr(attr,valattr)
    if is_verbose:
        print("{} = {}".format(attr,valattr))


# transcribe dimensions
for jdim in tb4com_dims.keys():
    if (tb4com_dims[jdim][0]=='dim'):
        dncdim = tb4com_dims[jdim][1]
        if(dncdim in dnc.dimensions.keys()):
            tbnc.createDimension(jdim,dnc.dimensions[dncdim].size)
        else:
            print("ERROR: Cannot find dimension {} in the input file.".format(dncdim))
    # later add alternative dimension creating options
    if  (tb4com_dims[jdim][0]=='fun'):
        dncdim = tb4com_dims[jdim][2][0]
        if(dncdim in dnc.dimensions.keys()):
            newval = tb_fun[tb4com_dims[jdim][1]](jdim,dnc.dimensions[dncdim].size)
            tbnc.createDimension(jdim,newval)
        else:
            print("ERROR: Cannot find dimension {} in the input file.".format(dncdim))


# transcribe variables
# dales tb order: ('time',<otherdimension>) 
# -> go by the list
for jvar in tb4com_vars.keys(): # tbvars.keys():
    tbnc.createVariable(jvar, tbvars[jvar]['dtype'], tbvars[jvar]['dimensions'], fill_value=def_fill_value)
    # fill depending on input type
    if ( tb4com_vars[jvar][0]=='var'): 
        # variable from variable
        # copy attributes
        for attr in dnc[tb4com_vars[jvar][1]].ncattrs():
            tbnc[jvar].setncattr(attr,dnc[tb4com_vars[jvar][1]].getncattr(attr))
        # check the dimensions of input
        dimsin = dnc[tb4com_vars[jvar][1]].dimensions
        if (dimsin[comsetup ['repeatdim'][0]]==comsetup ['repeatdim'][1]):
            # get the initial profile
            arrin = dnc[tb4com_vars[jvar][1]][:]
            # repeat it so it fits the dimensions
            arrout = np.squeeze( np.repeat(arrin, dnc.dimensions[comsetup ['repeatdim'][2]].size, axis = 0))
            # and insert values to output file
            tbnc[jvar][:]=arrout
        else: 
            # just remove singleton dimensions and transfer the array 
            tbnc[jvar][:]=np.squeeze(dnc[tb4com_vars[jvar][1]][:])
    elif ( tb4com_vars[jvar][0]=='att'):  
        # set default attributes
        for attr in atts_tbvars:
            tbnc[jvar].setncattr(attr,tbvars[jvar][attr])
        # variables from string attribute
        valin =dnc.getncattr(tb4com_vars[jvar][1])
        # convert to array of numerical values
        arrout = float(valin.split(' ')[0])*np.ones(tbnc[jvar].shape)
        # and insert values to the output file
        tbnc[jvar][:]=arrout
    elif ( tb4com_vars[jvar][0]=='fun'):
        # set default attributes
        for attr in atts_tbvars:
            tbnc[jvar].setncattr(attr,tbvars[jvar][attr])
            # later copy additional attributes?
        # does the function requires input ?
        if ( len(tb4com_vars[jvar][2])>0):
            #: later could be modified for more input arguments but not needed
            # get input variable to the function
            arrin = dnc[tb4com_vars[jvar][2][0]][:]
            # run the function
            arrout = tb_fun[tb4com_vars[jvar][1]](jvar,np.squeeze(arrin))
        else:
            arrout = tb_fun[tb4com_vars[jvar][1]](jvar)
        # and insert it
        tbnc[jvar][:]=arrout
    else:
        print("ERROR: The conversion option {} for varible {} not recognised {} = {}".format(tb4com_vars[jvar][0]),jvar)


# close files
dnc.close()
tbnc.close()


# get altitudes for the grid
dnc = Dataset(dephy_filename,'r')
zprof = dnc[prof4_com][:]

# write the grid file
write_prof_inp (zprof,fileout=out_dir+out_profname)








 





