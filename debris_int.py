import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import os
import string
import mcint
import random
import cld_fnc as cf
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=72, Om0=0.25)
from scipy.integrate import quad, dblquad, tplquad, nquad
from astropy.cosmology import z_at_value
import astropy.units as u
import time
from streamteam.util import get_pool

##############################################################
#general variables
jstokmkpcpersec = 3.241e-23
kmperkpc = 3.08567758*10.**16.
mperkpc  = 3.08567758*10.**19.
sminkg = 1.988435*10.**30.
secperyr = 3.15569e7
gee = 6.674*10**-11. #SI
gphys2 = gee/mperkpc*sminkg # in m^2 kpc / (M_sun s^2)
gphys3 = gee/(mperkpc**3.)*sminkg # in kpc^3/(M_sun s^2)

h_h0 = 0.72
rho_crit = 9.74*10**-27 * (h_h0/0.72)**2 # in kg/m^3

#halo parameters
mvir=1.77e12 #msun
rvir=389. #kpc
rhalo=24.6 #kpc
cc=rvir/rhalo
G=1
#rhalo = a
mhalo=mvir/(np.log(cc+1.)-cc/(cc+1.))

#mu choices
mu_shell = 1.
#debris_min_angle = np.pi/6.
debris_min_angle = 0.

ximin = 1e-3
ximax = 0.1
zmin = 0.1
zmax = 1.

domainsize = np.pi/2.*(ximax-ximin)*(zmax-zmin)

###############################################################

def sampler():
    while True:
        xi      = random.uniform(ximin,ximax)
        z       = random.uniform(zmin, zmax )
        vr      = random.uniform(0.00001, np.sqrt(2.)*.99999)
        vtheta  = random.uniform(0.00001, np.sqrt(2.-vr*vr))
        yield xi, z, vr, vtheta

#def intlimits(M, xi, z, v_r, v_theta):
	#	#integration limits
	#	ximin = 1e-5
	#	ximax = 0.1
	#	zmin = 0.1
	#	zmax = 1.
	#	rvir = (M*sminkg/(200.*rho_crit*(4./3.)*np.pi))**(1./3.)/mperkpc # virial radius in kpc
	#	mvir = m_encl_nfw(M,z,rvir)
	#	vcircvir = np.sqrt(gee*mvir*sminkg/rvir/mperkpc)
	#	vesc = vcircvir*np.sqrt(2)
	#	vrmin = 0.
	#	vtmin = 0.
	#	vrmax = vesc
	#	vtmax = 
	#

def integrand_2(x, args, verbose=False):
	M = args[0]
	df= args[1]

	xi 			= x[0] 
	z 			= x[1] 
	v_r 		= x[2] 
	v_theta 	= x[3] 

	dnmdxidz = 0.0104*(M/10**12.)**0.133 * xi**-1.995 * np.exp((xi/0.00972)**0.263) * (1+z)**0.0993
	rvir = (M*sminkg/(200.*rho_crit*(4./3.)*np.pi))**(1./3.)/mperkpc # virial radius in kpc
	c = halo_concentration(M, z)
	vcircvir = np.sqrt(gee*M*sminkg/rvir/mperkpc) #m/s
	#vescvir  = vcircvir*np.sqrt(2)
	e =  0.5*((v_r*vcircvir)**2.+(v_theta*vcircvir)**2) + cf.nfw_pot(M, z, rvir)
	l = v_theta*vcircvir*rvir*mperkpc

	#if verbose: print xi, z, v_r, v_theta, e, l

	msat = M*xi
	tinteract = cosmo.age(0).value - cosmo.age(z).value
	eps = 0.0000001
	local_l_start = np.array([[1.-eps,1.,1+eps],[1.-eps,1.,1+eps],[1.-eps,1.,1+eps]])
	local_e_start = np.array([[1.+eps,1.+eps,1.+eps],[1.,1.,1.],[1.-eps,1.-eps,1.-eps]])
	local_dps = np.zeros((3,3))
	local_trs = np.zeros((3,3))
	local_l = local_l_start*l
	local_e = local_e_start*e
	dlstep=(np.gradient(local_l)[1][1,1])
	destep=(np.gradient(local_e)[0][1,1])

	for iter1 in np.arange(3):
		for iter2 in np.arange(3):
			if ((iter1 ==1) | (iter2 == 1)):
				#time1 = time.time()
				rcirc                  = get_rcirc_e_nfw(         M,z,local_e[iter1,iter2], local_l[iter1,iter2]) 
				#print rcirc
				#print Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], rcirc
				apo, peri              = calc_apo_peri_bisect_nfw(M,z,local_e[iter1,iter2], local_l[iter1,iter2], rcirc=rcirc)
				#time3 = time.time()


				if iter1==1: local_dps[iter1,iter2] = delta_psi_nfw(           M,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
				if iter2==1: local_trs[iter1,iter2] = radial_period_nfw(       M,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
				if ((iter1 ==1) & (iter2 == 1)):
					r_t_onehalf        = r_t_onehalf_nfw(         M,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
					alpha              = delta_psi_nfw(           M,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=r_t_onehalf)
					
					#if verbose:print 'here', Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2]
					#if verbose:print 'rcirc', rcirc, 'apo', apo,'peri', peri
				#time4 = time.time()
				#time2 = time.time()

	ddpdl = np.gradient(local_dps,destep,dlstep)[1][1,1]
	dtrde = np.gradient(local_trs,destep,dlstep)[0][1,1]
	s, e_s, l_s = cf.scales_nfw(M, z, msat, l, peri)
	t_r = local_trs[1,1]
	dpsi= local_dps[1,1]
	omega_e = e_s*dtrde/t_r*dpsi*(tinteract/t_r)
	omega_e1 = np.min([alpha,omega_e])
	omega_l = l_s*ddpdl*(tinteract/t_r)
	mu = omega_l/(np.min([alpha,omega_e]))

	if df.lower() == 'uniform' : veldf = vel_func_uniform(v_r, v_theta) 
	if df.lower() == 'gaussian': veldf = vel_func_gaussian(v_r, v_theta)
	if df.lower() == 'benson':   veldf = vel_func_benson(v_r, v_theta, args[2])


	inte = dnmdxidz*veldf*(np.max([omega_e1,omega_l]) > debris_min_angle)
	inte_sh = inte * (mu > mu_shell)
	inte_st = inte * (mu < mu_shell)

	if verbose: print dnmdxidz, veldf, (mu > mu_shell), (np.max([omega_e1,omega_l]) > debris_min_angle), inte_sh
	

	#print time2-time1
	#print time4-time3
	return inte_sh, inte_st

def checkdf(df,*args,**keywords):
	plt.clf()
	x = np.linspace(0,2,200)
	y = np.linspace(0,2,200)
	X,Y = np.meshgrid(x,y)
	plt.subplot(aspect='equal')
	plt.imshow(df(X,Y, (keywords)), origin='lower', extent=[0,2,0,2])
	plt.plot(x, np.sqrt(2-x**2),c='k', linewidth=3)
	plt.colorbar()
	if df==vel_func_benson: 
		print 'integral in bound region', dblquad(df, 0.,np.sqrt(2), lambda x:0, lambda x:np.sqrt(2.-x**2), args = (keywords['redshift'],))
	elif df==vel_func_gaussian:
		print 'integral in bound region',   dblquad(df, 0.,np.sqrt(2), lambda x:0, lambda x:np.sqrt(2.-x**2), args = (keywords,))
	else: print 'integral in bound region', dblquad(df, 0.,np.sqrt(2), lambda x:0, lambda x:np.sqrt(2.-x**2))

	
def vel_func_benson(v_r, v_theta, redshift):
	if redshift=='0.0':
		a1 = 3.90
		a2 = 2.49
		a3 = 10.2
		a4 = .684
		a5 = .354
		a6 = 1.08
		a7 = .510
		a8 = .206
		a9 = .315

	if redshift=='0.5':
		a1 = 4.46
		a2 = 2.98
		a3 = 11.0
		a4 = 1.11
		a5 = .494
		a6 = 1.16
		a7 = .261
		a8 =-.279
		a9 = .331 

	if redshift=='1.0':
		a1 = 6.38
		a2 = 2.30
		a3 = 18.8
		a4 = .506
		a5 = -.0934
		a6 = 1.05
		a7 = .267
		a8 =-.154
		a9 = .157 

	b1 = a3*np.exp(-a4*(v_theta-a5)**2)
	b2 = a6*np.exp(-a7*(v_theta-a8)**2)
	val = a1*v_theta*np.exp(-a2*(v_theta-a9)**2-b1*(v_r-b2)**2)

	return val

def vel_func_uniform(v_r, v_theta):

	return v_r*0.+v_theta*0.+1./(np.pi*2./4.)

def vel_func_gaussian(vr, vtheta, *args):
	theta    = args[0].setdefault('theta',0.)
	center_r = args[0].setdefault('center_r', 1.)
	center_theta = args[0].setdefault('center_theta', 1.)
	sigma_r = args[0].setdefault('sigma_r', 1./3.)
	sigma_theta = args[0].setdefault('sigma_theta', 1./3.)

	#input vr and vtheta as fractions of the virial velocity
	norm = 1./(sigma_r*np.sqrt(2.*np.pi)) * 1./(sigma_theta*np.sqrt(2.*np.pi))

	a = np.cos(theta)**2/2/sigma_r**2 + np.sin(theta)**2/2/sigma_theta**2
	b = -np.sin(2*theta)/4/sigma_r**2 + np.sin(2*theta)/4/sigma_theta**2
	c = np.sin(theta)**2/2/sigma_r**2 + np.cos(theta)**2/2/sigma_theta**2
 
	df = norm * np.exp(-(a*(vr-center_r)**2. + 2.*b*(vr-center_r)*(vtheta-center_theta) + c*(vtheta-center_theta)**2.))

	return df

def halo_concentration(Mvir,z):
	a = 0.537+(1.025-0.537)*np.exp(-0.718*z**1.08)
	b = -0.097+0.024*z
	c = 10.**(a+b*np.log10(Mvir/10.**12.*h_h0))
	return c

def nshells(args, nsamples=10):
	time1 = time.time()
	nshell, nshell_err, nstream, nstream_err = mcint.integrate_2(integrand_2, sampler(), args = args, measure=domainsize, n=nsamples)
	time2 = time.time()
	return args, nshell, nshell_err, nstream, nstream_err

def dfreq():
	pool = get_pool(threads=4)
	ms = (10**((np.linspace(0,4,4)/8.+10)))
	args = [(i,'uniform') for i in ms]
	args = args + [(i,'gaussian') for i in ms]
	time1 = time.time()
	output = pool.map(nshells,args)
	time2 = time.time()
	pool.close()
	return output,time2-time1



if 0:



	clf()
	for i in np.arange(len(output)):
	    plt.errorbar(ms[i], output[i][0], yerr=output[i][1], c='r')
	    plt.errorbar(ms[i], output[i][2], yerr=output[i][3], c='b')

	clf()
	plt.errorbar(np.log10(output[output[:,1]=='gaussian',0].astype(float)), output[output[:,1]=='gaussian',2].astype(float), yerr=output[output[:,1]=='gaussian',3].astype(float),ls='-',  c='r', label='shell, gaussian')
	plt.errorbar(np.log10(output[output[:,1]=='uniform',0].astype(float)),  output[output[:,1]=='uniform',2].astype(float),  yerr=output[output[:,1]=='uniform',3].astype(float), ls='--',  c='r', label='shell, uniform')

	plt.errorbar(np.log10(output[output[:,1]=='gaussian',0].astype(float)), output[output[:,1]=='gaussian',4].astype(float), yerr=output[output[:,1]=='gaussian',5].astype(float),ls='-', c='b', label='stream, gaussian')
	plt.errorbar(np.log10(output[output[:,1]=='uniform',0].astype(float)),  output[output[:,1]=='uniform',4].astype(float),  yerr=output[output[:,1]=='uniform',5].astype(float), ls='--', c='b', label='stream, uniform')
	plt.legend(['shell, gaussian', 'shell, uniform', 'stream, gaussian', 'stream, uniform'], loc='upper left')
	plt.xlabel('$\mathrm{log \ M_{host}\ [M_\odot]}$')
	plt.ylabel('N [1]')

	clf()
	plt.errorbar(np.log10(output[output[:,1]=='uniform',0].astype(float)),   output[output[:,1]=='uniform', 2].astype(float)+output[output[:,1]=='uniform', 4].astype(float))
	plt.errorbar(np.log10(output[output[:,1]=='uniform',0].astype(float)),  (output[output[:,1]=='gaussian',2].astype(float)+output[output[:,1]=='gaussian',4].astype(float)))


def worker(i):
	x = 3**2.4 *6
	return x, x/3.


################lower accuracy integrators / root finders

if 1:	
	def apo_peri_func_nfw(r,M,z,e,l):
		n = (1./(r*mperkpc))**2. + 2*(cf.nfw_pot(M,z,r) - e) /(l)**2.
		return n 

	def calc_apo_peri_bisect_nfw(M, z, e, l, rcirc=25.):	
		#rcirc = get_rcirc_e(e,l)
		peri = scipy.optimize.bisect(apo_peri_func_nfw, .0001, rcirc, args = (M,z,e,l), maxiter = 20000)
		apo  = scipy.optimize.bisect(apo_peri_func_nfw, rcirc, 10000 , args = (M,z,e,l), maxiter = 20000)
		return apo, peri

	def dpsi_func_nfw(u, M, z, e, l):
		#takes peri and apo in m
		out = -1./(np.sqrt( 2.*( e - cf.nfw_pot(M, z, 1./(u*mperkpc))) - (l*u)**2. ))
		return out

	def delta_psi_nfw(M, z, e, l, apo = -1, peri= -1):
		#takes e and L in J and J*s, returns the radial period in Gyr
		if (apo < 0. and peri < 0): apo, peri = calc_apo_peri_bisect_nfw(e,l)
		part_d_psi, err = scipy.integrate.quad(dpsi_func_nfw, 1./(peri*mperkpc), 1./(apo*mperkpc), args = (M,z,e,l))
		return 2*l*part_d_psi

	def tr_func_nfw(r, M, z, e, l):
		#takes peri and apo in m
		out = 1./np.sqrt( 2.*( e - cf.nfw_pot(M,z,r/mperkpc) ) - ( l/(r) )**2. )
		return out

	def radial_period_nfw(M, z, e, l, apo = -1, peri= -1):
		#takes e and L in J and J*s, returns the radial period in Gyr
		if (apo < 0. and peri < 0): apo, peri = calc_apo_peri_bisect_nfw(M, z, e, l)
		half_t_r, err = scipy.integrate.quad(tr_func_nfw, peri*mperkpc, apo*mperkpc, args = (M, z, e, l))
		return (2.*half_t_r)/(secperyr * 10**9.)

	def r_t_onehalf_func_nfw(r,M,z,e,l,apo,t_r):
		thist = radial_period_nfw(M,z,e,l,apo=apo,peri=r)
		return thist - t_r/2.

	def r_t_onehalf_nfw(M, z, e, l, apo = -1, peri= -1):
		if (apo < 0. and peri < 0): apo, peri = calc_apo_peri_bisect_nfw(M, z, e, l)
		t_r = radial_period_nfw(M,z,e,l, apo=apo, peri=peri)
		r_t_onehalf = scipy.optimize.bisect(r_t_onehalf_func_nfw, peri, apo, args = (M,z,e,l,apo,t_r), maxiter = 20000)
		return r_t_onehalf

	def get_rcirc_e_func_nfw(r,M, z, ecirc, lcirc, verbose=False):
		if verbose:print  ecirc - (0.5*(gphys2*cf.m_encl_nfw(M,z,r)/r) + cf.nfw_pot(M,z,r))
		return ecirc - (0.5*(gphys2*cf.m_encl_nfw(M,z,r)/r) + cf.nfw_pot(M,z,r))

	def get_rcirc_e_nfw(M,z,ecirc,lcirc):
		rcirc = scipy.optimize.bisect(get_rcirc_e_func_nfw, 0.001, 5000., args = (M,z,ecirc,lcirc), maxiter=10000)
		return rcirc



#1.1343871651672042, 0.03844629554567381
#1.1372720818712685, 0.03741648527490802
#1.1201744632059341, 0.03677412106745073

