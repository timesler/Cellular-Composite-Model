/* Created by Language version: 6.2.0 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define _threadargscomma_ /**/
#define _threadargs_ /**/
 
#define _threadargsprotocomma_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gnabar _p[0]
#define gkbar _p[1]
#define gabar _p[2]
#define gcabar _p[3]
#define gkcbar _p[4]
#define idrk _p[5]
#define iak _p[6]
#define icak _p[7]
#define m_inf _p[8]
#define h_inf _p[9]
#define n_inf _p[10]
#define p_inf _p[11]
#define q_inf _p[12]
#define c_inf _p[13]
#define tau_m _p[14]
#define tau_h _p[15]
#define tau_n _p[16]
#define tau_p _p[17]
#define tau_q _p[18]
#define tau_c _p[19]
#define m_exp _p[20]
#define h_exp _p[21]
#define n_exp _p[22]
#define p_exp _p[23]
#define q_exp _p[24]
#define c_exp _p[25]
#define m _p[26]
#define h _p[27]
#define n _p[28]
#define p _p[29]
#define q _p[30]
#define c _p[31]
#define ena _p[32]
#define ek _p[33]
#define eca _p[34]
#define cao _p[35]
#define cai _p[36]
#define Dm _p[37]
#define Dh _p[38]
#define Dn _p[39]
#define Dp _p[40]
#define Dq _p[41]
#define Dc _p[42]
#define ina _p[43]
#define ik _p[44]
#define ica _p[45]
#define _g _p[46]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
#define _ion_cai	*_ppvar[6]._pval
#define _ion_eca	*_ppvar[7]._pval
#define _ion_cao	*_ppvar[8]._pval
#define _ion_ica	*_ppvar[9]._pval
#define _ion_dicadv	*_ppvar[10]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_evaluate_fct(void);
 static void _hoc_states(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_spike", _hoc_setdata,
 "evaluate_fct_spike", _hoc_evaluate_fct,
 "states_spike", _hoc_states,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gnabar_spike", "mho/cm2",
 "gkbar_spike", "mho/cm2",
 "gabar_spike", "mho/cm2",
 "gcabar_spike", "mho/cm2",
 "gkcbar_spike", "mho/cm2",
 "idrk_spike", "mA/cm2",
 "iak_spike", "mA/cm2",
 "icak_spike", "mA/cm2",
 0,0
};
 static double c0 = 0;
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 static double n0 = 0;
 static double p0 = 0;
 static double q0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"spike",
 "gnabar_spike",
 "gkbar_spike",
 "gabar_spike",
 "gcabar_spike",
 "gkcbar_spike",
 0,
 "idrk_spike",
 "iak_spike",
 "icak_spike",
 "m_inf_spike",
 "h_inf_spike",
 "n_inf_spike",
 "p_inf_spike",
 "q_inf_spike",
 "c_inf_spike",
 "tau_m_spike",
 "tau_h_spike",
 "tau_n_spike",
 "tau_p_spike",
 "tau_q_spike",
 "tau_c_spike",
 "m_exp_spike",
 "h_exp_spike",
 "n_exp_spike",
 "p_exp_spike",
 "q_exp_spike",
 "c_exp_spike",
 0,
 "m_spike",
 "h_spike",
 "n_spike",
 "p_spike",
 "q_spike",
 "c_spike",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 47, _prop);
 	/*initialize range parameters*/
 	gnabar = 0.04;
 	gkbar = 0.012;
 	gabar = 0.036;
 	gcabar = 0.002;
 	gkcbar = 5e-005;
 	_prop->param = _p;
 	_prop->param_size = 47;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 11, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 1);
 	_ppvar[6]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[7]._pval = &prop_ion->param[0]; /* eca */
 	_ppvar[8]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[9]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[10]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _spike_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	ion_reg("ca", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_prop_size(_mechtype, 47, 11);
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 spike C:/Users/tesler/Dropbox Data/PhD/MATLAB/Cellular Composite Model/Simulations/Spatial Shaping/Data/spike.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "HH style channels for spiking retinal ganglion cells";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int evaluate_fct(double);
static int states();
 
static int  states (  ) {
   evaluate_fct ( _threadargscomma_ v ) ;
   m = m + m_exp * ( m_inf - m ) ;
   h = h + h_exp * ( h_inf - h ) ;
   n = n + n_exp * ( n_inf - n ) ;
   p = p + p_exp * ( p_inf - p ) ;
   q = q + q_exp * ( q_inf - q ) ;
   c = c + c_exp * ( c_inf - c ) ;
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 
static int  evaluate_fct (  double _lv ) {
   double _la , _lb ;
 _la = ( - 0.6 * ( _lv + 30.0 ) ) / ( ( exp ( - 0.1 * ( _lv + 30.0 ) ) ) - 1.0 ) ;
   _lb = 20.0 * ( exp ( ( - 1.0 * ( _lv + 55.0 ) ) / 18.0 ) ) ;
   tau_m = 1.0 / ( _la + _lb ) ;
   m_inf = _la * tau_m ;
   _la = 0.4 * ( exp ( ( - 1.0 * ( _lv + 50.0 ) ) / 20.0 ) ) ;
   _lb = 6.0 / ( 1.0 + exp ( - 0.1 * ( _lv + 20.0 ) ) ) ;
   tau_h = 1.0 / ( _la + _lb ) ;
   h_inf = _la * tau_h ;
   _la = ( - 0.02 * ( _lv + 40.0 ) ) / ( ( exp ( - 0.1 * ( _lv + 40.0 ) ) ) - 1.0 ) ;
   _lb = 0.4 * ( exp ( ( - 1.0 * ( _lv + 50.0 ) ) / 80.0 ) ) ;
   tau_n = 1.0 / ( _la + _lb ) ;
   n_inf = _la * tau_n ;
   _la = ( - 0.006 * ( _lv + 90.0 ) ) / ( ( exp ( - 0.1 * ( _lv + 90.0 ) ) ) - 1.0 ) ;
   _lb = 0.1 * ( exp ( ( - 1.0 * ( _lv + 30.0 ) ) / 10.0 ) ) ;
   tau_p = 1.0 / ( _la + _lb ) ;
   p_inf = _la * tau_p ;
   _la = 0.04 * ( exp ( ( - 1.0 * ( _lv + 70.0 ) ) / 20.0 ) ) ;
   _lb = 0.6 / ( 1.0 + exp ( - 0.1 * ( _lv + 40.0 ) ) ) ;
   tau_q = 1.0 / ( _la + _lb ) ;
   q_inf = _la * tau_q ;
   _la = ( - 0.3 * ( _lv + 13.0 ) ) / ( ( exp ( - 0.1 * ( _lv + 13.0 ) ) ) - 1.0 ) ;
   _lb = 10.0 * ( exp ( ( - 1.0 * ( _lv + 38.0 ) ) / 18.0 ) ) ;
   tau_c = 1.0 / ( _la + _lb ) ;
   c_inf = _la * tau_c ;
   m_exp = 1.0 - exp ( - dt / tau_m ) ;
   h_exp = 1.0 - exp ( - dt / tau_h ) ;
   n_exp = 1.0 - exp ( - dt / tau_n ) ;
   p_exp = 1.0 - exp ( - dt / tau_p ) ;
   q_exp = 1.0 - exp ( - dt / tau_q ) ;
   c_exp = 1.0 - exp ( - dt / tau_c ) ;
    return 0; }
 
static void _hoc_evaluate_fct(void) {
  double _r;
   _r = 1.;
 evaluate_fct (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("spike", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 6, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 7, 0);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 8, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 9, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 10, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  c = c0;
  h = h0;
  m = m0;
  n = n0;
  p = p0;
  q = q0;
 {
   m = 0.0345 ;
   h = 0.8594 ;
   n = 0.1213 ;
   p = 0.0862 ;
   q = 0.2534 ;
   c = 0.0038 ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
  ek = _ion_ek;
  cai = _ion_cai;
  eca = _ion_eca;
  cao = _ion_cao;
 initmodel();
   }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ina = gnabar * m * m * m * h * ( v - ena ) ;
   idrk = gkbar * n * n * n * n * ( v - ek ) ;
   iak = gabar * p * p * p * q * ( v - ek ) ;
   icak = gkcbar * ( ( cai / 0.001 ) / ( 1.0 + ( cai / 0.001 ) ) ) * ( v - ek ) ;
   ik = idrk + iak + icak ;
   ica = gcabar * c * c * c * ( v - eca ) ;
   }
 _current += ina;
 _current += ik;
 _current += ica;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
  ek = _ion_ek;
  cai = _ion_cai;
  eca = _ion_eca;
  cao = _ion_cao;
 _g = _nrn_current(_v + .001);
 	{ double _dica;
 double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
  _dica = ica;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
  _ion_ik += ik ;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
 double _break, _save;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _break = t + .5*dt; _save = t;
 v=_v;
{
  ena = _ion_ena;
  ek = _ion_ek;
  cai = _ion_cai;
  eca = _ion_eca;
  cao = _ion_cao;
 { {
 for (; t < _break; t += dt) {
 error =  states();
 if(error){fprintf(stderr,"at line 84 in file spike.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 }   }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}
