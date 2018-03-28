from Timba.TDL import *
from Meow import Context

import random
import math
import pickle

beams = pickle.load(open('askap36_boresight_beams_simulated.p','rb'))

DEG = math.pi/180.0
ARCMIN = DEG/60
ARCSEC = DEG/3600


def askap_beam(E,lm,p):
  """
  Computes a beam for the given direction.
  'E' is output node.
  'lm' is direction (2-vector node).
  'p' is the current station.
  """
  ns = E.Subscope();

  ref_freq = 1.409e9
  ref_fwhm = 1.0159

  #ants = ['2','4','5','10', '12', '13', '14', '16', '24', '27', '28', '30']
  ants = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35']

  beam_key_xx = ants[int(p)]+','+str(beam_num+1)+',XX'
  beam_key_yy = ants[int(p)]+','+str(beam_num+1)+',YY'

  beam_params_xx = beams[beam_key_xx]
  beam_params_yy = beams[beam_key_yy]

  # beam parameters are all FWHM, whereas Gaussian equation
  # takes dispersions, hence the 2.355 divisions

  fscale = 3.0+((-2.0/ref_freq)*Meq.Freq())

  x0_x = beam_params_xx[0] * DEG
  y0_x = beam_params_xx[1] * DEG
  if use_ideal:
    ns.maj_x << (ref_fwhm/2.355) * fscale
    ns.min_x << (ref_fwhm/2.355) * fscale
    PA_x = 0.0
    A_x = 1.0
  else:
    ns.maj_x << (beam_params_xx[2]/2.355) * fscale
    ns.min_x << (beam_params_xx[3]/2.355) * fscale
    PA_x = beam_params_xx[4]
    A_x = beam_params_xx[5]

  x0_y = beam_params_yy[0] * DEG
  y0_y = beam_params_yy[1] * DEG
  if use_ideal:
    ns.maj_y << (ref_fwhm/2.355) * fscale
    ns.min_y << (ref_fwhm/2.355) * fscale
    PA_y = 0.0
    A_y = 1.0
  else:
    ns.maj_y << (beam_params_yy[2]/2.355) * fscale
    ns.min_y << (beam_params_yy[3]/2.355) * fscale
    PA_y = beam_params_yy[4]
    A_y = beam_params_yy[5]

  ns.lm0_x << Meq.Composer(x0_x,y0_x)
  ns.lm_x << (lm - ns.lm0_x) / DEG
  ns.l_x  << Meq.Selector(ns.lm_x, index=0)
  ns.m_x  << Meq.Selector(ns.lm_x, index=1)

  ns.lm0_y << Meq.Composer(x0_y,y0_y)
  ns.lm_y << (lm - ns.lm0_y) / DEG
  ns.l_y  << Meq.Selector(ns.lm_y, index=0)
  ns.m_y  << Meq.Selector(ns.lm_y, index=1)

  a_x = (Meq.Pow(ns.l_x,2))*((((math.cos(PA_x*DEG))**2.0)/(2.0*(ns.maj_x**2.0))) + (((math.sin(PA_x*DEG))**2.0)/(2.0*(ns.min_x**2.0))))
  b_x = 2.0*Meq.Multiply(ns.l_x,ns.m_x)*(((math.sin(2.0*PA_x*DEG))/(4.0*(ns.maj_x**2.0))) - ((math.sin(2.0*PA_x*DEG))/(4.0*(ns.min_x**2.0))))
  c_x = (Meq.Pow(ns.m_x,2))*((((math.sin(PA_x*DEG))**2.0)/(2.0*(ns.maj_x**2.0))) + (((math.cos(PA_x*DEG))**2.0)/(2.0*(ns.min_x**2.0))))

  a_y = (Meq.Pow(ns.l_y,2))*((((math.cos(PA_y*DEG))**2.0)/(2.0*(ns.maj_y**2.0))) + (((math.sin(PA_y*DEG))**2.0)/(2.0*(ns.min_y**2.0))))
  b_y = 2.0*Meq.Multiply(ns.l_y,ns.m_y)*(((math.sin(2.0*PA_y*DEG))/(4.0*(ns.maj_y**2.0))) - ((math.sin(2.0*PA_y*DEG))/(4.0*(ns.min_y**2.0))))
  c_y = (Meq.Pow(ns.m_y,2))*((((math.sin(PA_y*DEG))**2.0)/(2.0*(ns.maj_y**2.0))) + (((math.cos(PA_y*DEG))**2.0)/(2.0*(ns.min_y**2.0))))

  ns.E_x << A_x * Meq.Exp(-1.0*(a_x+b_x+c_x))
  ns.E_y << A_y * Meq.Exp(-1.0*(a_y+b_y+c_y))


  E << Meq.Matrix22(ns.E_x, 0, 0, ns.E_y);
  return E;



def compute_jones (Jones,sources,stations=None,pointing_offsets=None,**kw):
  """Computes beam gain for a list of sources.
  The output node, will be qualified with either a source only, or a source/station pair
  """;
  stations = stations or Context.array.stations();
  ns = Jones.Subscope();
  # are pointing errors configured?
  if pointing_offsets:
    # create nodes to compute actual pointing per source, per antenna
    for p in stations:
      for src in sources:
        lm = ns.lm(src.direction,p) << src.direction.lm() + pointing_offsets(p);
        beam_model(Jones(src,p),lm,p);
  # no pointing errors
  else:
      for src in sources:
        for p in stations:
          beam_model(Jones(src,p),src.direction.lm(),p);
  return Jones;


# def compute_jones_tensor (Jones,sources,stations=None,lmn=None,
#                           label="beam",pointing_offsets=None,inspectors=[],
#                           **kw):
#   """Computes beam gain tensor for a list of sources.
#   The output node, will be qualified with either a source only, or a source/station pair
#   """;
#   stations = stations or Context.array.stations();
#   ns = Jones.Subscope();
  
#   # figure out beam model
#   for beam_model in MODELS:
#     if globals().get('use_%s'%beam_model.label,None):
#       break;
#   else:
#     raise RuntimeError,"no beam model selected";
  
#   if getattr(beam_model,'compute_tensor',None) is None:
#     return None;
    
#   beam_model.prepare(ns);
  
#   # see if sources have a "beam_lm" or "_lm_ncp" attribute
#   lmsrc = [ src.get_attr("beam_lm",None) or src.get_attr("_lm_ncp",None) for src in sources ];
  
#   # if all source have the attribute, create lmn tensor node (and override the lmn argument)
#   if all([lm is not None for lm in lmsrc]):
#     lmn = ns.lmnT << Meq.Constant(lmsrc);
    
#   # if lmn tensor is not set for us, create a composer
#   if lmn is None:
#     lmn = ns.lmnT << Meq.Composer(dims=[0],*[ src.direction.lm() for src in sources ]);

#   # create station tensors
#   for ip,p in enumerate(stations):
#     beam_model.compute_tensor(Jones(p),lmn,pointing_offsets and pointing_offsets(p),ip);
          
#   return Jones;

_model_option = TDLCompileOption('beam_model',"Beam model",
  [askap_beam]
);

_askap_option_menu = TDLCompileMenu('Select ASKAP beam',
  TDLOption('beam_num', "Beam number",
    [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35],
    more=int,doc="ASKAP beam to simulate"),
  TDLMenu('Use ideal beam',toggle="use_ideal",doc="Click this to use circular Gaussians with FWHM of 1 deg."))

def _show_option_menus (model):
  _askap_option_menu.show(model==askap_beam);

_model_option.when_changed(_show_option_menus);
