import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.integrate
import os
import string
import mcint
import random
from scipy.integrate import quad, dblquad, tplquad, nquad


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

#halo parameters
mvir=1.77e12 #msun
rvir=389. #kpc
rhalo=24.6 #kpc
cc=rvir/rhalo
G=1
#rhalo = a
mhalo=mvir/(np.log(cc+1.)-cc/(cc+1.))

##############################################################
#utils

def printclear(somestring, end=False):
	from sys import stdout
	stdout.write("\r" + somestring)
	stdout.flush()
	if end: 
		stdout.write("\n")
	return


##############################################################
#functions to read snap/com files

def read_snap(dir, snap):
	#read snap files and put them into physical units
	pars = np.loadtxt(dir + 'SCFPAR', dtype='string', usecols=(0,1))
	ru = float(pars[16,0]) #kpc
	mu = float(pars[17,0]) #msun

	gphys1 = 4.49*10**(-6) # in kpc^3/(M_sun Gyr^2)
	tu = np.sqrt(ru**3./(gphys1 * mu)) # in Gyr
	vu = np.sqrt(mu*gphys1/ru) * 3.086*10.**16. * 10.**(-9.) * 3.171*10**-8# in (kpc/Gyr ->) km/s
	eu = vu**2. * 10.**6. # in J per unit mass

	if snap < 10: snapf = 'SNAP00' + str(snap)
	if snap >= 10 and snap < 100: snapf = 'SNAP0' + str(snap)
	if snap >= 100: snapf = 'SNAP' + str(snap)

	tsnap = float(np.loadtxt(dir+snapf, usecols=(0,1),dtype='string')[0,1])
	data = np.loadtxt(dir+snapf, skiprows = 1)
	# mass x y z vx vy vz ep potext tub

	if data.shape[1] == 10:
		print 'Remnant still bound'
		m = data[:,0] * mu
		x = data[:,1] * ru
		y = data[:,2] * ru
		z = data[:,3] * ru
		vx = data[:,4] * vu
		vy = data[:,5] * vu
		vz = data[:,6] * vu
		ep = data[:,7] * eu
		potex = data[:,8] * eu
		tub = data[:,9] * tu
		return m, x, y, z, vx, vy, vz, ep, potex, tub, tsnap*tu

	if data.shape[1] == 8:
		print 'Remnant disrupted'
		m = data[:,0] * mu
		x = data[:,1] * ru
		y = data[:,2] * ru
		z = data[:,3] * ru
		vx = data[:,4] * vu
		vy = data[:,5] * vu
		vz = data[:,6] * vu
		potex = data[:,7] * eu
		return m, x, y, z, vx, vy, vz, 0.*x, potex, 0.*x, tsnap*tu


def read_snap_sg(dir, snap):
	#read snap files and put them into physical units
	pars = np.loadtxt(dir + 'SCFPAR', dtype='string')
	ru = float(pars[16,0]) #kpc
	mu = float(pars[17,0]) #msun

	gphys1 = 4.49*10**(-6) # in kpc^3/(M_sun Gyr^2)
	tu = np.sqrt(ru**3./(gphys1 * mu)) # in Gyr
	vu = np.sqrt(mu*gphys1/ru) * 3.086*10.**16. * 10.**(-9.) * 3.171*10**-8# in (kpc/Gyr ->) km/s
	eu = vu**2. * 10.**6. # in J per unit mass

	if snap < 10: snapf = 'SNAP00' + str(snap)
	if snap >= 10 and snap < 100: snapf = 'SNAP0' + str(snap)
	if snap >= 100: snapf = 'SNAP' + str(snap)

	data = np.loadtxt(dir+snapf, skiprows = 1)
	# mass x y z vx vy vz ep potext tub

	m = data[:,0] * mu
	x = data[:,1] * ru
	y = data[:,2] * ru
	z = data[:,3] * ru
	vx = data[:,4] * vu
	vy = data[:,5] * vu
	vz = data[:,6] * vu
	ep = data[:,7] * eu

	return m, x, y, z, vx, vy, vz, ep

def read_tsnap_sg(dir, snap):
	#read snap files and put them into physical units
	pars = np.loadtxt(dir + 'SCFPAR', dtype='string')
	ru = float(pars[16,0]) #kpc
	mu = float(pars[17,0]) #msun

	gphys1 = 4.49*10**(-6) # in kpc^3/(M_sun Gyr^2)
	tu = np.sqrt(ru**3./(gphys1 * mu)) # in Gyr
	vu = np.sqrt(mu*gphys1/ru) * 3.086*10.**16. * 10.**(-9.) * 3.171*10**-8# in (kpc/Gyr ->) km/s
	eu = vu**2. * 10.**6. # in J per unit mass

	if snap < 10: snapf = 'TSNAP00' + str(snap)
	if snap >= 10 and snap < 100: snapf = 'TSNAP0' + str(snap)
	if snap >= 100: snapf = 'TSNAP' + str(snap)

	data = np.loadtxt(dir+snapf)
	# mass x y z vx vy vz ep potext tub

	m = data[:,0] * mu
	x = data[:,1] * ru
	y = data[:,2] * ru
	z = data[:,3] * ru
	vx = data[:,4] * vu
	vy = data[:,5] * vu
	vz = data[:,6] * vu
	ep = data[:,7] * eu

	return m, x, y, z, vx, vy, vz, ep


def read_com(dir):
	#read COM file and give data physical units
	pars = np.loadtxt(dir + 'SCFPAR', dtype='string', usecols=(0,1))
	ru = float(pars[16,0]) #kpc
	mu = float(pars[17,0]) #msun

	gphys1 = 4.49*10**(-6) # in kpc^3/(M_sun Gyr^2)
	tu = np.sqrt(ru**3./(gphys1 * mu)) # in Gyr
	vu = np.sqrt(mu*gphys1/ru) * kmperkpc * 10.**(-9.) * 3.171*10**-8. # in kpc/Gyr -> km/s

	com = np.loadtxt(dir + 'SCFCEN')
	#time dt x y z vx vy vz

	t = com[:,0] * tu
	dt = com[:,1] * tu
	x_c = com[:,2] * ru
	y_c = com[:,3] * ru
	z_c = com[:,4] * ru
	vx_c = com[:,5] * vu
	vy_c = com[:,6] * vu
	vz_c = com[:,7] * vu

	return t, dt, x_c, y_c, z_c, vx_c, vy_c, vz_c

def com_init(dir):
	initial_vals = np.loadtxt(dir+'SCFPAR', dtype='string', skiprows = 18)
	i_spt = [float(initial_vals[0,0]), float(initial_vals[0,1]), float(initial_vals[0,2])]
	i_vel = [float(initial_vals[1,0]), float(initial_vals[1,1]), float(initial_vals[1,2])]

	return i_spt[0], i_spt[1], i_spt[2], i_vel[0], i_vel[1], i_vel[2]



##############################################################
#derived quantities

def get_Ls(x,y,z,vx,vy,vz):
	#x in kpc, vx in km/s, result is J*s per unit mass
	L_x = (y*vz - z*vy)*mperkpc*1000.#*3.086*10.**22.        #yz - zy
	L_y = (z*vx - x*vz)*mperkpc*1000.#*3.086*10.**22.        #zx - xz
	L_z = (x*vy - y*vx)*mperkpc*1000.#*3.086*10.**22.        #xy - yx
	L_mag = np.sqrt(L_x**2. + L_y**2. + L_z**2.)

	return L_x, L_y, L_z, L_mag

def apo_snaps(dir):
	t, dt, x_c, y_c, z_c, vx_c, vy_c, vz_c = read_com(dir)
	return 0 

def radial(x,y,z,vx,vy,vz):
	r = np.sqrt(x**2. + y**2. + z**2.)
	vr = (x*vx + y*vy + z*vz)/abs(r)
	theta = np.arccos(z/r)
	phi = np.arctan2(y,x)
	b = 90. - theta*180./np.pi
	l =  phi*180./np.pi

	return r, theta, phi, vr, l, b

def radial_s(x,y,z,vx,vy,vz):
	x = x+8.33
	r = np.sqrt(x**2. + y**2. + z**2.)
	vr = (x*vx + y*vy + z*vz)/abs(r)
	theta = np.arccos(z/r)
	phi = np.arctan2(y,x)
	b = 90. - theta*180./np.pi
	l =  phi*180./np.pi

	return r, theta, phi, vr, l, b



##############################################################
#general functions of the potential

def halo_pot(r):
	#r in kpc
	p = r/rhalo

	gphys2 = gee/mperkpc*sminkg # in m^2 kpc / (M_sun s^2)
	
	phi0=gphys2*mhalo/rhalo
	
	phi=-(phi0/p)*np.log(p+1.) # in m^2/s^2

	return phi

def dphi_dr(r):
	p=r/rhalo
	dphi_dr_ = mhalo*gee*sminkg/(mperkpc**2.) * ((np.log(1.+p)/r**2.) - 1./(rhalo*r*(1+p)))
	return dphi_dr_

def dsqrphi_drsqr(r):
	#mhalo = 4 pi rho0 a^3
	#rho0 = mhalo/(4.*np.pi*rhalo**3.)
	p=r/rhalo
	dsqrphi_drsquare = mhalo*gee*sminkg/(mperkpc**3.) * ((-2.*np.log(1.+p)/r**3.) + (2.*rhalo+3.*r)/(rhalo**2.*r**2.*(1+p)**2.))
	#return in J/m^2/kg
	return dsqrphi_drsquare


def rho(r):
	#return the density at r
	rho0 = mhalo/(4.*np.pi*rhalo**3.)
	p=r/rhalo
	rho = rho0 / (p * (1+p)**2.)
	return rho

def mass_shell(r):
	#return the mass of a shell at r, for calculating enclosed mass (this matches m_encl)
	return rho(r)*r**2.*4.*np.pi

def m_encl(r):
	#returns the halo mass enclosed by r (in kpc) in solar masses
	p=r/rhalo
	m_encl=mhalo*(np.log(p+1.)-p/(p+1.))
	return m_encl

def scales(m_sat, l_sat, r_scale):
	#returns s = r_j/r, l_s = l/l_c, and e_s = gm/r
	#m_sat in msun, r_scale in kpc
	s = (m_sat/(3.*m_encl(r_scale)))**(1./3.)
	e_s = 2. * s * gphys2*m_encl(r_scale)/r_scale
	l_s = 3.732 * s*l_sat
	return s, e_s, l_s

def king_scales(m_sat, l_sat, peri): # e_sat):
	#returns s = r_j/r, l_s = l/l_c, and e_s = gm/r
	#m_sat in msun, r_scale in kpc
	w_peri = l_sat/(peri*mperkpc)**2.
	r_lim = (gee*sminkg*m_sat/(w_peri**2 - dsqrphi_drsqr(peri)))**(1./3.)
	s = r_lim/mperkpc/peri
	e_s = 2.*s * peri * mperkpc *dphi_dr(peri)
	l_s = 3.732*s*l_sat
	return s, e_s, l_s


##############################################################
#orbit properties as f(l,e)

def get_vs_func(vs, ra, rp):
	return (vs*1000 - np.sqrt(2*(halo_pot(ra)-halo_pot(rp)) + (rp/ra*vs*1000.)**2.))


def get_vs(ra, rp):
	vp = scipy.optimize.bisect(get_vs_func,1000.,1., maxiter = 10000, args=(ra,rp))
	va = vp*rp/ra
	return va, vp

def get_rcirc_l_func(r, ecirc, lcirc):
	return lcirc - r*mperkpc*np.sqrt(gphys2*m_encl(r)/r)

def get_rcirc_l(ecirc,lcirc):
	rcirc = scipy.optimize.bisect(get_rcirc_l_func, 0.01,500., args = (ecirc,lcirc), maxiter=10000, tol = 1e-9)
	return rcirc

def get_rcirc_e_func(r,ecirc, lcirc):
	return ecirc - (0.5*(np.sqrt(gphys2*m_encl(r)/r))**2. + halo_pot(r))

def get_rcirc_e(ecirc,lcirc):
	rcirc = scipy.optimize.bisect(get_rcirc_e_func, 0.01, 500., args = (ecirc,lcirc), maxiter=10000)
	return rcirc


def apo_peri_func(r,e,l):
	n = (1./(r*mperkpc))**2. + 2*(halo_pot(r) - e) /(l)**2.
	return n 


def calc_apo_peri(e, l):
	#apo/peri ~= ps[0]*ls**4 + ps[1]*ls**3 + ps[2]*ls**2 + ps[3]*ls + ps[4], ls = l/lcirc
	peri = scipy.optimize.newton(apo_peri_func, 1.,  args = (e,l), maxiter = 10000, tol=1e-6)
	apo = scipy.optimize.newton(apo_peri_func, 50., args = (e,l), maxiter = 10000, tol=1e-6)
	if (apo-peri)/apo < .05: 
		print 'warning 1: apo and peri possibly the same:', apo, peri
		peri = scipy.optimize.newton(apo_peri_func, 1.,  args = (e,l), maxiter = 10000, tol=1e-6)
		apo = scipy.optimize.newton(apo_peri_func, rvir*2., args = (e,l), maxiter = 10000, tol=1e-6)
		if (apo-peri)/apo < .05: 
			print 'warning 2: apo and peri possibly the same:', apo, peri
			ps = [147.07104192, -450.17504672,  501.41350007, -248.32220742, 50.72046458]
			qs = [ 0.51896956,  0.25131472,  0.03008319]
			rcirc = get_rcirc_e(e,l)
			lcirc = rcirc*mperkpc*np.sqrt(gphys2*m_encl(rcirc)/rcirc)
			ls = l/lcirc
			#apo/peri ~= ps[0]*ls**4 + ps[1]*ls**3 + ps[2]*ls**2 + ps[3]*ls + ps[4], ls = l/lcirc
			peri_guess = rcirc*(qs[0]*ls**2 + qs[1]*ls + qs[2])
			apo_guess = peri_guess*(ps[0]*ls**4 + ps[1]*ls**3 + ps[2]*ls**2 + ps[3]*ls + ps[4])
			peri = scipy.optimize.newton(apo_peri_func, peri_guess*.5,  args = (e,l), maxiter = 10000, tol=1e-6)
			apo = scipy.optimize.newton(apo_peri_func, apo_guess*1.25, args = (e,l), maxiter = 10000, tol=1e-6)
			
			print apo, peri
	return apo, peri

def calc_apo_peri_bisect(e, l, rcirc=25.):	
	#rcirc = get_rcirc_e(e,l)
	peri = scipy.optimize.bisect(apo_peri_func, .0001, rcirc, args = (e,l), maxiter = 20000)
	apo  = scipy.optimize.bisect(apo_peri_func, rcirc, 1000,  args = (e,l), maxiter = 20000)
	return apo, peri


def tr_func(r, e, l):
	#takes peri and apo in m
	out = 1./np.sqrt( 2.*( e - halo_pot(r/mperkpc) ) - ( l/(r) )**2. )
	return out

def radial_period(e, l, apo = -1, peri= -1):
	#takes e and L in J and J*s, returns the radial period in Gyr
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri(e,l)
	half_t_r, err = scipy.integrate.quad(tr_func, peri*mperkpc, apo*mperkpc, args = (e,l))
	return (2.*half_t_r)/(secperyr * 10**9.)

def r_t_onehalf_func(r,e,l,apo,t_r):
	thist = radial_period(e,l,apo=apo,peri=r)
	return thist - t_r/2.

def r_t_onehalf_fidu(e, l, apo = -1, peri= -1):
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri_bisect(e, l)
	t_r = radial_period(e,l, apo=apo, peri=peri)
	r_t_onehalf = scipy.optimize.bisect(r_t_onehalf_func, peri, apo, args = (e,l,apo,t_r), maxiter = 20000)
	return r_t_onehalf

###########################################################
#actions

def j_r_func(r, e, l):
	#takes peri and apo in m
	out = np.sqrt( 2.*( e - halo_pot(r/mperkpc) ) - ( l/(r) )**2. )
	return out

def j_r(e, l, apo = -1, peri= -1):
	#takes e and L in J and J*s, returns the radial period in Gyr
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri_bisect(e,l)
	integr, err = scipy.integrate.quad(j_r_func, peri*mperkpc, apo*mperkpc, args = (e,l))
	return (1./np.pi)*integr

def calc_H_func(h, jr, l):
	return jr - j_r(h,l)


def calc_H(jr, l, start = -2.8e10):
	H =scipy.optimize.newton(calc_H_func, start, args = (jr, l), maxiter = 3000, tol = 1e-9)
	return H



###########################################################
#angular functions, substituting u = 1/r to remove singularities

def dpsi_func(u, e, l):
	#takes peri and apo in m
	out = -1./(np.sqrt( 2.*( e - halo_pot(1./(u*mperkpc))) - (l*u)**2. ))
	return out


def delta_psi(e, l, apo = -1, peri= -1, rcirc=25.):
	#takes e and L in J and J*s, returns the radial period in Gyr
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri_bisect(e,l,rcirc=rcirc)
	part_d_psi, err = scipy.integrate.quad(dpsi_func, 1./(peri*mperkpc), 1./(apo*mperkpc), args = (e,l))
	return 2*l*part_d_psi


def azi_period(e, l, apo = -1, peri= -1, rcirc=25.):
	#returns t_psi in Gyr
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri_bisect(e,l,rcirc=rcirc)
	t_psi = 2.*np.pi*radial_period(e,l, apo=apo, peri=peri)/(abs(delta_psi(e,l, apo=apo, peri=peri)))
	return t_psi

def omega_precess(e, l, apo=-1, peri=-1):
	# in rad per Gyr
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri(e,l)
	omega_p = (delta_psi(e,l, apo=apo, peri=peri) - 2*np.pi)/radial_period(e,l,apo=apo, peri=peri)
	return omega_p


def simple_omegas_emperical(e, l, r_sat, v_sat, rcirc=25.):
	#given e and l, return vperi/rperi, vapo/rapo, vcirc/rcirc
	apo, peri = calc_apo_peri(e,l)
	vcirc = np.sqrt(gphys2*m_encl(rcirc)/rcirc)/1000. #km/s
	omega_circ = vcirc/(rcirc*kmperkpc)
	omega_apo  = min(v_sat)/max(r_sat*kmperkpc)
	omega_peri = max(v_sat)/min(r_sat*kmperkpc)

	return omega_circ, omega_apo, omega_peri # in rad/s



#######################################################
#######################################################
#find orbital resonances using E since the radial period is dominated by it

def radial_res_func(e, l, tr_sat, res_period):
	out = res_period - radial_period(e,l)
	print e, out
	return out

def radial_resonance(e, l, m):
	tr_sat = radial_period(e,l)
	#divde by m instead of multiply since omega goes like 1/t_r
	res_period = tr_sat/m
	res_energy = scipy.optimize.newton(radial_res_func, e, args = (l, tr_sat, res_period))
	return res_energy

def azi_res_func(e, l, t_sat, res_period):
	out = res_period - azi_period(e,l)
	print e, out
	return out

def azi_resonance(e, l, m):
	t_sat = azi_period(e,l)
	#divde by m instead of multiply since omega goes like 1/t_r
	res_period = t_sat/m
	res_energy = scipy.optimize.newton(azi_res_func, e, args = (l, t_sat, res_period))
	return res_energy


#######################################################
#leapfrog integrator
kvgee=6.67e-8
msun=1.989e33
cmpkpc=3.085678e21
secperyr=3.1536e7


def int_init(x,y,z,vx,vy,vz,dt = 0.02, reduc=False,	ru = .136, mu = 2.5e6):

	tu=np.sqrt((cmpkpc*ru)**3/(msun*mu*kvgee))
	vu = (ru*cmpkpc*1.e-5)/tu

	if reduc:

		x=x/ru
		y=y/ru
		z=z/ru
		vx=vx/vu
		vy=vy/vu
		vz=vz/vu

	return x,y,z,vx,vy,vz,ru,mu,tu,vu,dt



def accp_tri(x,y,z,vx,vy,vz, model = 'lm10', dt = 0.02,ru=0.136, mu=2.5e6):

	#C -----------------------------------------------------------------------
	#C       Three Component Galaxy model:
		#C
		#C               Miyamoto-Nagai Disk + Hernquist Spheroid + Logarithmic Halo
		#C
		#C
		#C       Features:
		#C
		#C               Nearly flat rotation curve from 1 - 30 kpc
		#C               Disk z scale height of 200 pc
		#C               kappa_z has correct radial dependance
		#C
		#C
		#C       Mass inside solar circle: (roughly)  
		#C
		#C
		#C               Disk 42%
		#C               Spheroid 28%
		#C               Halo 30%
		#C --------------------------------------------------------------
		#C Convert parameters to appropriate units....

	x,y,z,vx,vy,vz,ru,mu,tu,vu,dt = int_init(x,y,z,vx,vy,vz, ru=ru, mu=mu)

	a=6.5/ru
	b=.26/ru
	c=12./ru
	rs=.7/ru
	vcirc2=220.**2/vu/vu
	GM = 8.887*vcirc2/ru
	GMs = 3.0*vcirc2/ru
	vh2 = 121.858*121.858/vu/vu

	#if model=='newberg':
	#	#spherical halo from LJM 05
	#	vh2 = 114.*114./vu/vu
	#	q1=1.


	if model.lower() =='lm10':
		#triaxial halo from LM 10
		q1=1.38
		q2=1.0
		qz=1.36
		#phi in radians (97 degrees)
		phi=1.692969
		C1 = (np.cos(phi)/q1)**2+(np.sin(phi)/q2)**2
		C2 = (np.cos(phi)/q2)**2+(np.sin(phi)/q1)**2
		C3 = 2.*np.sin(phi)*np.cos(phi)*(1./q1/q1-1./q2/q2)

	if model.lower() =='lm10_spherical':
		#triaxial halo from LM 10
		q1=1.0
		q2=1.0
		qz=1.0
		#phi in radians (97 degrees)
		phi=1.692969
		C1 = (np.cos(phi)/q1)**2+(np.sin(phi)/q2)**2
		C2 = (np.cos(phi)/q2)**2+(np.sin(phi)/q1)**2
		C3 = 2.*np.sin(phi)*np.cos(phi)*(1./q1/q1-1./q2/q2)

	if model.lower() =='lm10_prolate':
		#triaxial halo from LM 10
		q1=1.0
		q2=1.0
		qz=1.36
		#phi in radians (97 degrees)
		phi=1.692969
		C1 = (np.cos(phi)/q1)**2+(np.sin(phi)/q2)**2
		C2 = (np.cos(phi)/q2)**2+(np.sin(phi)/q1)**2
		C3 = 2.*np.sin(phi)*np.cos(phi)*(1./q1/q1-1./q2/q2)

	if model.lower() =='lm10_oblate':
		#triaxial halo from LM 10
		q1=1.0
		q2=1.0
		qz=0.68
		#phi in radians (97 degrees)
		phi=1.692969
		C1 = (np.cos(phi)/q1)**2+(np.sin(phi)/q2)**2
		C2 = (np.cos(phi)/q2)**2+(np.sin(phi)/q1)**2
		C3 = 2.*np.sin(phi)*np.cos(phi)*(1./q1/q1-1./q2/q2)

	if model.lower() =='chaos':
		#triaxial halo from LM 10
		q1=1.242
		q2=1.0
		qz=1.367
		#phi in radians (97 degrees)
		phi=1.692969
		C1 = (np.cos(phi)/q1)**2+(np.sin(phi)/q2)**2
		C2 = (np.cos(phi)/q2)**2+(np.sin(phi)/q1)**2
		C3 = 2.*np.sin(phi)*np.cos(phi)*(1./q1/q1-1./q2/q2)

	#compute acceleration due to spheroid

	r2 = x*x+y*y
	z2 = z*z
	rad = np.sqrt(r2+z2)
	tsrad = GMs/(rad+rs)**2/rad
	phis = -GMs/(rad+rs)

	ax = -tsrad*x
	ay = -tsrad*y
	az = -tsrad*z

	
	#compute acceleration due to disk
	
	sqz2b2 = np.sqrt(z2 + b*b)
	tdr = GM/(r2 + (a + sqz2b2)**2)**1.5
	tdz = tdr*(a/sqz2b2 + 1.)
	phim = -GM/np.sqrt(r2+(a+sqz2b2)**2)

	ax = ax -tdr*x                         
	ay = ay -tdr*y                          
	az = az -tdz*z


	#compute acceleration due to halo

	#if (q1 == 1.):
	#	thrad = 2.*vh2/(r2 + z2 + c*c)
	#	phih = vh2*np.sqrt(r2+z2+c*c)
	#	ax = ax -thrad*x
	#	ay = ay -thrad*y
	#	az = az -thrad*z
	#else:
	xx=x
	yy=y
	zz=z

	re2=C1*xx*xx+C2*yy*yy+C3*xx*yy+zz*zz/(qz*qz)
	phih=vh2*np.log(re2+c*c)

	fac=vh2/(re2+c*c)
	dpdx=fac*(2.*C1*xx+C3*yy)
	dpdy=fac*(2.*C2*yy+C3*xx)
	dpdz=fac*(2.*zz/qz/qz)

	ax = ax - dpdx
	ay = ay - dpdy
	az = az - dpdz


	return ax,ay,az


def stepsys(x,y,z,vx,vy,vz,ax,ay,az,dt):
	x,y,z    = cf.steppos(x,y,z,vx,vy,vz,dt)
	ax,ay,az = cf.get_accel(x,y,z,vx,vy,vz)
	vx,vy,vz = cf.stepvel(vx,vy,vz,ax,ay,az,dt)

def start_vel(vx,vy,vz,ax,ay,az,dt):
	vx = vx+0.5*ax*dt
	vy = vy+0.5*ay*dt
	vz = vz+0.5*az*dt
	return vx,vy,vz

def stop_vel(vx,vy,vz,ax,ay,az,dt):
	vx = vx-0.5*ax*dt
	vy = vy-0.5*ay*dt
	vz = vz-0.5*az*dt
	return vx,vy,vz

def stepvel(vx,vy,vz,ax,ay,az,dt):
	vx = vx+ax*dt
	vy = vy+ay*dt
	vz = vz+az*dt
	return vx,vy,vz

def steppos(x,y,z,vx,vy,vz,dt):
	x = x+vx*dt
	y = y+vy*dt
	z = z+vz*dt
	return x,y,z

def get_accel(x,y,z,vx,vy,vz, mu=2.5e8, ru=0.2, vu=73.3):
	#halo parameters
	mvir=1.77e12 #msun
	rvir=389. #kpc
	rhalo=24.6 #kpc
	cc=rvir/rhalo
	G=1.
	gee=6.67e-8
	ggg=gee*msun/cmpkpc/1.e5/1.e5

	#rhalo = a
	mhalo=mvir/(np.log(cc+1.)-cc/(cc+1.))
	phi0=ggg*mhalo/rhalo
	mhalo=mhalo/mu
	rhalo=rhalo/ru
	phi0=phi0/vu/vu

	r=np.sqrt(x**2.+y**2.+z**2.)
	p=r/rhalo
	mr=mhalo*(np.log(p+1.)-p/(p+1.))
	aror=-(G*mr/r/r)/r
	phi=-(phi0/p)*np.log(p+1.)
	ax = aror*x
	ay = aror*y
	az = aror*z

	return ax,ay,az




def get_accel_s(x,y,z,vx,vy,vz):

	r = np.sqrt(x**2.+y**2.+z**2.)
	ax = -gee*m_encl(r/mperkpc)*sminkg*x/r**3.#halo_pot(r/mperkpc)*x/r**2
	ay = -gee*m_encl(r/mperkpc)*sminkg*y/r**3.#halo_pot(r/mperkpc)*y/r**2
	az = -gee*m_encl(r/mperkpc)*sminkg*z/r**3.#halo_pot(r/mperkpc)*z/r**2

	return ax,ay,az


#for scfpy debugging
def pb():
	global x,y,z,vx,vy,vz, ax,ay,az
	print 'x0 before:  ', x[0]
	print 'y0 before:  ', y[0]
	print 'z0 before:  ', z[0]
	print 'vx0 before: ', vx[0]
	print 'vy0 before: ', vy[0]
	print 'vz0 before: ', vz[0]
	print 'ax0 before: ', ax[0]
	print 'ay0 before: ', ay[0]
	print 'az0 before: ', az[0]

	return

def pa():
	global x,y,z,vx,vy,vz, ax,ay,az
	print 'x0 after:  ', x[0]
	print 'y0 after:  ', y[0]
	print 'z0 after:  ', z[0]
	print 'vx0 after: ', vx[0]
	print 'vy0 after: ', vy[0]
	print 'vz0 after: ', vz[0]
	print 'ax0 after: ', ax[0]
	print 'ay0 after: ', ay[0]
	print 'az0 after: ', az[0]

	return



##############################################################
#the particle class that will hold each snap's data

class particles:
	def init(self, dire, snap, do_jr=False, rel = True):
		q = snap+1
		self.m, self.x, self.y, self.z, self.vx, self.vy, self.vz, self.ep, self.potex, self.tub, self.tsnap = read_snap(dire, snap)
		self.r, self.theta, self.phi, self.vr, self.l, self.b = radial(self.x, self.y, self.z, self.vx, self.vy, self.vz)
		self.r_s, self.theta_s, self.phi_s, self.vr_s, self.l_s, self.b_s = radial_s(self.x, self.y, self.z, self.vx, self.vy, self.vz)
		self.l_x, self.l_y, self.l_z, self.l_mag = get_Ls(self.x, self.y, self.z, self.vx, self.vy, self.vz)
		self.te  = 0.5*1000.*1000.*((self.vx)**2.+(self.vy)**2+(self.vz)**2.) + self.potex
		self.nbods = len(self.x)

		if rel: 
			self.t, self.dt, self.x_c, self.y_c, self.z_c, self.vx_c, self.vy_c, self.vz_c = read_com(dire)
			self.l_x_c, self.l_y_c, self.l_z_c, self.l_mag_c = get_Ls( self.x_c, self.y_c, self.z_c, self.vx_c, self.vy_c, self.vz_c)
			self.e_c = 0.5*1000.*1000.*((self.vx_c[q])**2.+(self.vy_c[q])**2+(self.vz_c[q])**2.) + halo_pot(np.sqrt(self.x_c[q]**2. + self.y_c[q]**2.+ self.z_c[q]**2.))
			self.de = self.te - self.e_c
			self.dl = self.l_mag - self.l_mag_c[q]
			self.apo_c, self.peri_c = calc_apo_peri_bisect(self.e_c, self.l_mag_c[q])
			self.s, self.e_s, self.l_s = scales(np.sum(self.m), self.l_mag_c[q], self.peri_c)

		if do_jr:
			self.j_r = np.zeros(self.nbods)
			if rel: self.j_r_c = j_r(self.e_c, self.l_mag_c[q])
			for i in np.arange(self.nbods):
				try:
					self.j_r[i] = j_r(self.te[i], self.l_mag[i])
				except: 
					self.j_r[i] = -1
			print '%2.1f percent of particles have bad j_r'%(100.*np.sum(self.j_r<0)/self.nbods)




def lbr_to_xyz(l,b,r, sun_pos = -8):	
	from astropy.coordinates import GalacticCoordinates, ICRSCoordinates
	from astropy import units as u
	from astropy.coordinates import Distance

	ccoord = GalacticCoordinates(b=b, l=l, unit=(u.degree, u.degree), distance=Distance(r, u.kpc))

	x = (ccoord.x+sun_pos)
	y = ccoord.y
	z = ccoord.z

	return x,y,z



def yetissh():
	import paramiko
	import getpass

	ssh = paramiko.SSHClient()
	ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
	localpath = '/scratch/hendel/clouds/code/nfw/clouds/'
	path = "/vega/astro/users/dah2154/nfw/clouds/"
	server = 'yetisubmit.cc.columbia.edu'
	username = 'dah2154'
	password = getpass.getpass()
	ssh.connect(server, username=username, password=password)
	sftp = ssh.open_sftp()

	return ssh, sftp


def downloadsim(localpath,yetipath):
	#localpath = '/media/backup1/hendel/sims/clouds/rc45/'
	#yetipath  = '/vega/astro/users/dah2154/nfw/clouds/rc45/'
	ssh,sftp = yetissh()
	massdirs = ['M2.5e+06','M2.5e+07','M2.5e+08','M2.5e+09']
	ldirs = ['0.05','0.10','0.20','0.30','0.40','0.50','0.60','0.70','0.80','0.90','0.95','1.00']


	nsims = len(ldirs)*len(massdirs)
	counter = 0

	os.system('mkdir ' + localpath)

	for iter1 in massdirs:
		os.system('mkdir ' + localpath + iter1)

		for iter2 in ldirs:
			thispath = iter1 + '/L' + iter2 + '/'
			os.system('mkdir ' + localpath + thispath)

			for files in sftp.listdir(yetipath+thispath):
				sftp.get(yetipath+thispath+files, localpath+thispath+files)
			
			counter = counter+1
			print '%3.0f of %3.0f sims copied'%(counter,nsims)

	sftp.close()
	ssh.close()




def check_e_consv(dire):
	plt.plot(np.loadtxt(dire+'SCFLOG')[:,-1]*100.)
	plt.ylabel('total energy change [percent of ep0]')
	plt.xlabel('noutlog')
	return

def mass_loss(dire, snap=-1, nts= 200.):
	files = np.sort(os.listdir(dire))

	for f in files:
		if (string.find(f,"SNAP")!=-1): 
			lastsnap = f

	lastsnapn = string.atoi(string.split(lastsnap, 'SNAP')[1])
	data = particles()
	data.init(dire,lastsnapn, rel=False)
	ts = np.arange(nts)/nts*np.max(data.tub)
	ms = ts*0.
	for t in np.arange(200):
		ms[t] = np.sum(data.m[(data.tub==0) | (data.tub > ts[t])])

	plt.plot(ts,ms)

def showapos(dire):
	pars = np.loadtxt(dire + 'SCFPAR', dtype='string', usecols=(0,1))
	nlog = float(pars[3,0])
	nbod = float(pars[4,0])
	t, dt, x_c, y_c, z_c, vx_c, vy_c, vz_c = read_com(dire)
	r = np.sqrt(x_c**2+y_c**2+z_c**2)
	plt.plot(t,r)
	for i in np.arange(len(r)-2)+1:
		if ((r[i]>r[i-1]) & (r[i]>r[i+1])):
			plt.scatter(t[i],r[i], s=10, c = 'r')
			plt.text(t[i], r[i], '%2.2f'%((i+1)/(nbod/nlog)))




#########################################################################
#########################################################################
#Functions to calculate cloud frequency
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=72, Om0=0.25)
#get age at a redshift
#cosmo.age(1).value

def nfw_pot(Mvir, z, r, h=0.72):
	a = 0.537+(1.025-0.537)*np.exp(-0.718*z**1.08)
	b = -0.097+0.024*z
	c = 10.**(a+b*np.log10(Mvir/10.**12.*h))
	#print a,b,c
	#rho_crit = 3 H **2 / (8 pi G)
	rho_crit = 9.74*10**-27 * (h/0.72)**2 # in kg/m^3
	rvir = (Mvir*sminkg/(200.*rho_crit*(4./3.)*np.pi))**(1./3.)/mperkpc # virial radius in kpc
	rscale = rvir/c
	#print rho_crit,rvir,rscale
	#r in kpc
	p = r/rscale
	mhalo=Mvir/(np.log(c+1.)-c/(c+1.))
	phi0=gphys2*mhalo/rscale
	#print mhalo, phi0
	phi=-(phi0/p)*np.log(p+1.) # in m^2/s^2
	return phi

def apo_peri_func_nfw(r,M,z,e,l):
	n = (1./(r*mperkpc))**2. + 2*(nfw_pot(M,z,r) - e) /(l)**2.
	return n 

def calc_apo_peri_bisect_nfw(M, z, e, l, rcirc=25.):	
	#rcirc = get_rcirc_e(e,l)
	peri = scipy.optimize.bisect(apo_peri_func_nfw, .0001, rcirc, args = (M,z,e,l), maxiter = 20000)
	apo  = scipy.optimize.bisect(apo_peri_func_nfw, rcirc, 1000, args = (M,z,e,l), maxiter = 20000)
	return apo, peri

def dpsi_func_nfw(u, M, z, e, l):
	#takes peri and apo in m
	out = -1./(np.sqrt( 2.*( e - nfw_pot(M, z, 1./(u*mperkpc))) - (l*u)**2. ))
	return out

def delta_psi_nfw(M, z, e, l, apo = -1, peri= -1):
	#takes e and L in J and J*s, returns the radial period in Gyr
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri_bisect_nfw(e,l)
	part_d_psi, err = scipy.integrate.quad(dpsi_func_nfw, 1./(peri*mperkpc), 1./(apo*mperkpc), args = (M,z,e,l))
	return 2*l*part_d_psi

def tr_func_nfw(r, M, z, e, l):
	#takes peri and apo in m
	out = 1./np.sqrt( 2.*( e - nfw_pot(M,z, r/mperkpc) ) - ( l/(r) )**2. )
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

def m_encl_nfw(Mvir,z,r,h=0.72):
	#returns the halo mass enclosed by r (in kpc) in solar masses
	a = 0.537+(1.025-0.537)*np.exp(-0.718*z**1.08)
	b = -0.097+0.024*z
	c = 10.**(a+b*np.log10(Mvir/10.**12.*h))
	#print a,b,c
	#rho_crit = 3 H **2 / (8 pi G)
	rho_crit = 9.74*10**-27 * (h/0.72)**2 # in kg/m^3
	rvir = (Mvir*sminkg/(200.*rho_crit*(4./3.)*np.pi))**(1./3.)/mperkpc # virial radius in kpc
	rscale = rvir/c
	#print rho_crit,rvir,rscale
	#r in kpc
	mhalo=Mvir/(np.log(c+1.)-c/(c+1.))
	p=r/rscale
	m_encl=mhalo*(np.log(p+1.)-p/(p+1.))
	return m_encl

def get_rcirc_e_func_nfw(r,M, z, ecirc, lcirc, verbose=False):
	if verbose:print  ecirc - (0.5*(gphys2*m_encl_nfw(M,z,r)/r) + nfw_pot(M,z,r))
	return ecirc - (0.5*(gphys2*m_encl_nfw(M,z,r)/r) + nfw_pot(M,z,r))

def get_rcirc_e_nfw(M,z,ecirc,lcirc):
	rcirc = scipy.optimize.bisect(get_rcirc_e_func_nfw, 0.001, 5000., args = (M,z,ecirc,lcirc), maxiter=10000)
	return rcirc

def scales_nfw(M,z,m_sat, l_sat, r_scale):
	#returns s = r_j/r, l_s = l/l_c, and e_s = gm/r
	#m_sat in msun, r_scale in kpc
	mencl = m_encl_nfw(M,z,r_scale)
	s = (m_sat/(3.*mencl))**(1./3.)
	e_s = 2. * s * gphys2*mencl/r_scale
	l_s = 3.732 * s * l_sat
	return s, e_s, l_s

def orbitalparameters(Mgal,msat,rp,eta,z,specific=True, h=0.72):
	if specific: msat=1.
	#ms in msun, rp  & eta are unitless
	rho_crit = 9.74*10**-27 * (h/0.72)**2 # in kg/m^3
	rvir = (Mgal*sminkg/(200.*rho_crit*(4./3.)*np.pi))**(1./3.)/mperkpc # virial radius in kpc
	rp = rp*mperkpc*rvir #rp to m

	rm = Mgal*msat/(Mgal+msat) #reduced mass
	ecc = np.sqrt(1.-eta**2 ) #eccentricity
	print rp/mperkpc, rvir, ecc
	L = np.sqrt(    rp*(1+ecc)*gee*sminkg*msat*Mgal*rm       )/rm
	E = L**2./(2.*(rvir*mperkpc)**2.) - gee*sminkg*Mgal/(rvir*mperkpc) 
	#convert E from point mass -> nfw
	print L**2./(2.*rp**2.), - gee*sminkg*Mgal/(rp) , nfw_pot(Mgal,z,rp/mperkpc)
	
	#E = (ecc**2 -1)/(1+ecc)*(gee*sminkg*Mgal*msat/2./rp)
	#E = E + gee*sminkg*Mgal/(rvir) + nfw_pot(Mgal,z,rvir/mperkpc)
	return E, L


def omegaeomegal(Mgal,xi,z,e,l, verbose=False, fidu=False):
	msat = Mgal*xi
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

	if not fidu:
		for iter1 in np.arange(3):
			for iter2 in np.arange(3):
				if ((iter1 ==1) | (iter2 == 1)):
					
					rcirc                  = get_rcirc_e_nfw(         Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2]) 
					if verbose: print rcirc
					#print Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], rcirc
					apo, peri              = calc_apo_peri_bisect_nfw(Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], rcirc=rcirc)
					local_dps[iter1,iter2] = delta_psi_nfw(           Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
					local_trs[iter1,iter2] = radial_period_nfw(       Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
					if ((iter1 ==1) & (iter2 == 1)):
						r_t_onehalf = r_t_onehalf_nfw(Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], apo = apo, peri = peri)
						alpha = delta_psi_nfw(Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=r_t_onehalf)
						
						if verbose:print 'here', Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2]
						if verbose:print 'rcirc', rcirc, 'apo', apo,'peri', peri

		ddpdl = np.gradient(local_dps,destep,dlstep)[1][1,1]
		dtrde = np.gradient(local_trs,destep,dlstep)[0][1,1]
		s, e_s, l_s = scales_nfw(Mgal,z,msat, l, peri)
		t_r = local_trs[1,1]
		dpsi= local_dps[1,1]
		omega_e = e_s*dtrde/t_r*dpsi*(tinteract/t_r)
		omega_e1 = np.min([alpha,omega_e])
		omega_l = l_s*ddpdl*(tinteract/t_r)
		mu = omega_l/(np.min([alpha,omega_e]))

	if fidu:
		for iter1 in np.arange(3):
			for iter2 in np.arange(3):
				if ((iter1 ==1) | (iter2 == 1)):					
					rcirc                  = get_rcirc_e(         local_e[iter1,iter2], local_l[iter1,iter2]) 
					if verbose: print rcirc
					apo, peri              = calc_apo_peri_bisect(local_e[iter1,iter2], local_l[iter1,iter2], rcirc=rcirc)
					local_dps[iter1,iter2] = delta_psi(           local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
					local_trs[iter1,iter2] = radial_period(       local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
					if ((iter1 ==1) & (iter2 == 1)):
						#alpha = delta_psi(local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=rcirc)
						rt_onehalf = r_t_onehalf_fidu(local_e[iter1,iter2], local_l[iter1,iter2], apo = apo, peri = peri)
						alpha = delta_psi(local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=rt_onehalf)

						if verbose:print 'here', Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2]
						if verbose:print 'rcirc', rcirc, 'apo', apo,'peri', peri

		ddpdl = np.gradient(local_dps,destep,dlstep)[1][1,1]
		dtrde = np.gradient(local_trs,destep,dlstep)[0][1,1]
		s, e_s, l_s = scales(msat, l, peri)
		t_r = local_trs[1,1]
		dpsi= local_dps[1,1]
		omega_e = e_s*dtrde/t_r*dpsi*(tinteract/t_r)
		omega_l = l_s*ddpdl*(tinteract/t_r)
		omega_e1 = np.min([alpha,omega_e])
		mu = omega_l/(np.min([alpha,omega_e]))

	return omega_e, omega_e1, omega_l, alpha, mu

def classifyer(n):
	from astropy.cosmology import z_at_value
	import astropy.units as u

	d=particles()
	mus = np.zeros(n)
	omegaes = np.zeros(n)
	omegae1s = np.zeros(n)
	omegals = np.zeros(n)
	alphas = np.zeros(n)
	usedrcircs =np.zeros(n)
	usedmasses =np.zeros(n)
	usedls = np.zeros(n)
	tsnaps = np.zeros(n)
	classifications = np.zeros(n)

	for i in np.arange(n):
		rootdir = '/media/backup1/hendel/sims/clouds/'
		randrcirc = np.random.choice([25.,50.,100.])
		randmass = np.random.choice(['M2.5e+06','M2.5e+07','M2.5e+08'])#,'M2.5e+09'])
		randl    = np.random.choice(['0.05','0.10','0.20','0.30','0.40','0.50','0.60','0.70','0.80','0.90','0.95'])
		randt    = np.random.randint(10,101)

		vcirc = np.sqrt(gee*m_encl(randrcirc)*sminkg/(randrcirc*mperkpc))/1000. #km/s
		lcirc = randrcirc*mperkpc*vcirc*1000. #J*s
		ecirc = halo_pot(randrcirc) + 0.5*(vcirc*1000.)**2

		d.init(rootdir+'rc%i/'%(int(randrcirc))+randmass+'/L'+randl+'/',randt, rel=False)
		plt.clf()
		plt.subplot(1,1,1,aspect='equal')
		plt.scatter(d.x,d.y, alpha = 0.15, s=2, edgecolor='none', c='k')
		#plt.scatter(d.x[0:30000],d.y[0:30000], alpha = 0.2, s=2, edgecolor='none', c='k')

		classification = raw_input("Type 1 for stream, 2 for shell, 3 for intermediate, 4 for not debris, 0 for unknown:")
		classifications[i] = int(classification)
		omegaes[i], omegae1s[i], omegals[i], alphas[i], mus[i] = omegaeomegal(1.77e11,float(string.split(randmass,'M')[1])/1.77e11,z_at_value(cosmo.age,cosmo.age(0)-d.tsnap*u.Gyr),ecirc,lcirc*float(randl), verbose=False, fidu=True)
		usedrcircs[i] = randrcirc
		usedmasses[i] = float(string.split(randmass,'M')[1])
		usedls[i] = float(randl)
		tsnaps[i] = d.tsnap
		print float(i)/float(n)*100.

	
	dout = np.transpose(np.vstack([usedrcircs, usedmasses, usedls, tsnaps, omegaes, omegae1s, omegals, alphas, mus, classifications]))

	np.savetxt('classifyer_output5.txt', dout)
	return dout
	#return omegaes, omegae1s, omegals, alphas, mus, classifications
	
	#plt.hist(mus[classifications==1], bins=np.linspace(0,20,10))
	#plt.hist(mus[classifications==2], bins=np.linspace(0,20,10))
	#plt.hist(mus, bins=np.linspace(0,20,10), histtype='step')


def check_angles():
	from astropy.cosmology import z_at_value
	import astropy.units as u

	d=particles()
	rootdir = '/media/backup1/hendel/sims/clouds/'
	randrcirc = np.random.choice([25.,50.,100.])
	randmass = np.random.choice(['M2.5e+06','M2.5e+07','M2.5e+08','M2.5e+09'])
	randl    = np.random.choice(['0.05','0.10','0.20','0.30','0.40','0.50','0.60','0.70','0.80','0.90','0.95'])
	randt    = np.random.randint(1,101)

	vcirc = np.sqrt(gee*m_encl(randrcirc)*sminkg/(randrcirc*mperkpc))/1000. #km/s
	lcirc = randrcirc*mperkpc*vcirc*1000. #J*s
	ecirc = halo_pot(randrcirc) + 0.5*(vcirc*1000.)**2

	d.init(rootdir+'rc%i/'%(int(randrcirc))+randmass+'/L'+randl+'/',randt, rel=False)
	plt.clf()
	plt.subplot(1,1,1,aspect='equal')
	#plt.scatter(d.x,d.y, alpha = 0.1, s=2, edgecolor='none', c='k')
	
	rotation_angle = d.tsnap / azi_period(ecirc, lcirc*float(randl), rcirc = randrcirc) * 2.*np.pi
	#print azi_period(ecirc, lcirc*float(randl), rcirc = randrcirc), d.tsnap 
	#print rotation_angle
	omegaes, omegae1s, omegals, alphas, mus = omegaeomegal(1.77e11,float(string.split(randmass,'M')[1])/1.77e11,z_at_value(cosmo.age,cosmo.age(0)-d.tsnap*u.Gyr),ecirc,lcirc*float(randl), verbose=False, fidu=True)

	newx =   d.x*np.cos(rotation_angle) + d.y*np.sin(rotation_angle)
	newy = - d.x*np.sin(rotation_angle) + d.y*np.cos(rotation_angle)
	plt.scatter(newx,newy, alpha = 0.1, s=2, edgecolor='none', c='k')
	plt.scatter(d.x, d.y, alpha = 0.1, s=2, edgecolor='none', c='r')
	print omegaes, omegae1s, omegals, alphas, mus 

	plt.plot([0,50],[0,  50*np.tan(alphas/2.)], color = 'r', linewidth=3)
	plt.plot([0,50],[0, -50*np.tan(alphas/2.)], color = 'r', linewidth=3)

	plt.plot([0,50],[0,  50*np.tan(omegae1s/2.)], color = 'g', linewidth=2)
	plt.plot([0,50],[0, -50*np.tan(omegae1s/2.)], color = 'g', linewidth=2)

	plt.plot([0,50],[0,  50*np.tan(omegals/2.)], color = 'b', linewidth=2)
	plt.plot([0,50],[0, -50*np.tan(omegals/2.)], color = 'b', linewidth=2)


#old, using rp and eta (yuck)
	#def omegaeomegal(Mgal,xi,z,rp,eta):
	#	msat = Mgal*xi
	#	e,l = orbitalparameters(Mgal,msat,rp,eta,z)
	#	tinteract = cosmo.age(0).value - cosmo.age(z).value
	#	eps = 0.0000001
	#	local_l_start = np.array([[1.-eps,1.,1+eps],[1.-eps,1.,1+eps],[1.-eps,1.,1+eps]])
	#	local_e_start = np.array([[1.+eps,1.+eps,1.+eps],[1.,1.,1.],[1.-eps,1.-eps,1.-eps]])
	#	local_dps = np.zeros((3,3))
	#	local_trs = np.zeros((3,3))
	#	local_l = local_l_start*l
	#	local_e = local_e_start*e
	#	dlstep=(np.gradient(local_l)[1][1,1])
	#	destep=(np.gradient(local_e)[0][1,1])
	#
	#	for iter1 in np.arange(3):
	#		for iter2 in np.arange(3):
	#			if ((iter1 ==1) | (iter2 == 1)):
	#				
	#				rcirc                  = get_rcirc_e_nfw(         Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2]) 
	#				print rcirc
	#				apo, peri              = calc_apo_peri_bisect_nfw(Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], rcirc=rcirc)
	#				
	#				local_dps[iter1,iter2] = delta_psi_nfw(           Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
	#				local_trs[iter1,iter2] = radial_period_nfw(       Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
	#				if ((iter1 ==1) & (iter2 == 1)): 
	#					alpha = delta_psi_nfw(Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=rcirc)
	#					
	#					print 'here', Mgal,z,local_e[iter1,iter2], local_l[iter1,iter2]
	#					print 'rcirc', rcirc, 'apo', apo,'peri', peri
	#
	#	ddpdl = np.gradient(local_dps,destep,dlstep)[1][1,1]
	#	dtrde = np.gradient(local_trs,destep,dlstep)[0][1,1]
	#	s, e_s, l_s = scales(msat, l, peri)
	#	t_r = local_trs[1,1]
	#	dpsi= local_dps[1,1]
	#	omega_e = e_s*dtrde/t_r*dpsi*(tinteract/t_r)
	#	omega_l = l_s*ddpdl*(tinteract/t_r)
	#	mu = omega_l/(np.min([alpha,omega_e]))
	#	return omega_e, omega_l, alpha, mu


def gaussveldist(vr, vtheta, center_r = 1., center_theta = 1., sigma_r = 1./3., sigma_theta = 1./3.):
	#input vr and vtheta as fractions of the virial velocity
	norm = 1./(sigma_r*np.sqrt(2.*np.pi)) * 1./(sigma_theta*np.sqrt(2.*np.pi))
	df = norm * np.exp(-(vr-center_r)**2./(2.*sigma_r**2.)) * np.exp(-(vtheta-center_theta)**2./(2.*sigma_theta**2.))
	return df

def orbitalparameters_vs(M, z, vr, vtheta,h=0.72):
	#vr and vtheta as fractions of vcirc
	#rho_crit = 3 H **2 / (8 pi G)
	rho_crit = 9.74*10**-27 * (h/0.72)**2 # in kg/m^3
	rvir = (M*sminkg/(200.*rho_crit*(4./3.)*np.pi))**(1./3.)/mperkpc # virial radius in kpc
	mvir = m_encl_nfw(M,z,rvir)
	#print M, mvir
	vcircvir = np.sqrt(gee*mvir*sminkg/rvir/mperkpc)
	#print vcircvir, rvir
	vesc = vcircvir*np.sqrt(2)
	e = nfw_pot(M, z, rvir) + 0.5*(vr*vr+vtheta*vtheta)*vcircvir*vcircvir
	l = vtheta*vcircvir*rvir*mperkpc
	return e, l, vcircvir

def find_vesc_nfw(M, z, h=0.72):
	#vr and vtheta as fractions of vcirc
	#rho_crit = 3 H **2 / (8 pi G)
	rho_crit = 9.74*10**-27 * (h/0.72)**2 # in kg/m^3
	rvir = (M*sminkg/(200.*rho_crit*(4./3.)*np.pi))**(1./3.)/mperkpc # virial radius in kpc
	mvir = m_encl_nfw(M,z,rvir)
	vcircvir = np.sqrt(gee*mvir*sminkg/rvir/mperkpc)
	vesc = vcircvir*np.sqrt(2)
	return vesc


def dnmdxidz(M,xi,z):
	#M in Msun
	dnmdxidz = 0.0104*(M/10**12.)**0.133 * xi**-1.995 * np.exp((xi/0.00972)**0.263) * (1+z)**0.0993
	return dnmdxidz

def dnmdxidz_mlast(xi,z,M):
	#M in Msun
	dnmdxidz = 0.0104*(M/10**12.)**0.133 * xi**-1.995 * np.exp((xi/0.00972)**0.263) * (1+z)**0.0993
	return dnmdxidz

def shell_frac_integrand(vtheta,vr,z,xi, M, mu_shell=0.275):
	vesc = find_vesc_nfw(M,z)
	e, l, vcircvir= orbitalparameters_vs(M, z, vr/vesc, vtheta/vesc)
	omega_e, omega_e1, omega_l, alpha, mu = omegaeomegal(M,xi,z,e,l)
	print vtheta,vr,z,xi,M
	print  gaussveldist(vr,vtheta), dnmdxidz(M,xi,z), (mu > mu_shell), (np.max([omega_e1,omega_l]) > np.pi/2.)
	return gaussveldist(vr,vtheta)*dnmdxidz(M,xi,z)*(mu > mu_shell)*(np.max([omega_e1,omega_l]) > np.pi/2.)

def stream_frac_integrand(xi, z, vr, vtheta, M, mu_shell=0.275):
	e, l, vcircvir, vesc = orbitalparameters_vs(M, z, vr, vtheta)
	omega_e, omega_e1, omega_l, alpha, mu = omegaeomegal(M,xi,z,e,l)
	print  gaussveldist(vr,vtheta), dnmdxidz(M,xi,z), (mu < mu_shell), (np.max([omega_e1,omega_l]) > np.pi/2.)
	return gaussveldist(vr,vtheta)*dnmdxidz(M,xi,z)*(mu < mu_shell)*(np.max([omega_e1,omega_l]) > np.pi/2.)


def veldistfunc(M,ximin,ximax,zmin,zmax):
	return shellfraction

def mergerfunc(M,ximin,ximax,zmin,zmax):
	dnmdm = dblquad(dnmdxidz_mlast,zmin,zmax, lambda x: ximin, lambda x: ximax, args=(M,))
	return dnmdm

def ndebris(M, mu_shell=0.275):
	ximin = 1e-5
	ximax = 0.1
	zmin = 0.
	zmax = 1.
	def xilim(M):
		return [ximin,ximax]
	def zlim(xi,M):
		return [zmin, zmax ]
	def vrlim(z,xi,M):
		print 'vrlim', [0.00001, find_vesc_nfw(M,z)], z, M, xi
		return [0.00001, find_vesc_nfw(M,z)]
	def vthetalim(vr,z,xi,M):
		print M, xi, z, vr
		print 'vthetalim', [0.00001, np.sqrt(find_vesc_nfw(M,z)**2.-vr**2.)]
		return [0.00001, np.sqrt(find_vesc_nfw(M,z)**2.-vr**2.)]

	nshells = nquad(shell_frac_integrand, [vthetalim, vrlim, zlim, xilim], args=(M,))

	return nshells, nstreams


#def mcint_shells

#def dfs(M,z,rp,eta,h=0.72):
#	M_star = 10**(12.42 - 1.56*z +0.38*z**2.)/h
#	C_0 = 3.380*(1 + 0.567*(M/M_star)**0.152)
#	C_1 = 0.242*(1 + 2.360*(M/M_star)**0.108)
#	R_0 = 3.140*(1 + 0.152*((1+z)**-4 * M/M_star)**0.410)
#	R_1 = 0.450*(1 - 0.395*((1+z)**-4 * M/M_star)**0.109)
#	dfdeta = C_0*eta**1.05*(1-eta)**(C_1)
#	dfdrp = R_0*np.exp(-(rp/R_1)**0.85)
#	return dfdrp*dfdeta
#
#def shell_integrand(xi,z,rp,eta,Mgal,mu_shell):
#	omega_e, omega_l, alpha, mu = omegaeomegal(Mgal,xi,z,rp,eta)
#	return dfs(Mgal,z,rp,eta)*dnmdxidz(Mgal,xi,z)*(mu > mu_shell)*(np.max([omega_e,omega_l]) > np.pi/2.)
#
#def nshells(Mgal, ximin=1e-5, ximax=1e-1, zmin=0., zmax = 1., rpmin=0.0001, rpmax=1., etamin=0., etamax=1., mu_shell=0.4):
#	nshells = scipy.integrate.nquad(shell_integrand, [[zmin,zmax],[ximin,ximax],[rpmin, rpmax],[etamin,etamax]], args = (Mgal,mu_shell))
#	return nshells

##
	#def M_star(z,h=0.72):
	#	return 10**(12.42 - 1.56*z +0.38*z**2.)/h
	#
	#def C0(M,z):
	#	C_0 = 3.38*(1 + 0.567*(M/M_star(z))**0.152)
	#	return C_0
	#
	#def C1(M,z):
	#	C_1 = 0.242*(1 + 2.36*(M/M_star(z))**0.108)
	#	return C_1
	#
	#def R0(M,z):
	#	R_0 = 3.14*(1 + 0.152*((1+z)**-4 * M/M_star(z))**0.410)
	#	return R_0
	#
	#def R1(M,z):
	#	R_1 = 0.450*(1 - 0.395*((1+z)**-4 * M/M_star(z))**0.109)
	#	return R_1
	#
	#def dfdeta(M,z,eta, h=0.72):
	#	dfdeta = C0(M,z)*eta**1.05*(1-eta)**(C1(M,z))
	#	return dfdeta
	#
	#def dfdrp(M,z,rp):
	#	dfdrp = R0(M,z)*np.exp(-(rp/R1(M,z))**0.85)
	#	return dfdrp



#def dndpars(M,z,rp,eta,xi):
#	return dfdrp(M,z,rp)*dfdeta(M,z,eta)*dnmdxidz(M,xi,z)




###make video

#for i in np.arange(500)+1:
#	plt.clf()
#	data = cf.particles()
#	data.init('../7_26/10e6/0.700000/',i)
#	plt.scatter(data.x,data.y, alpha = 0.2, color = 'k', edgecolor='none', s=6)
#	plt.xlabel('x [kpc]')
#	plt.ylabel('y [kpc]')
#	plt.xlim([-40,40])
#	plt.ylim([-40,40])
#	plt.xlabel('$x [kpc]$')
#	plt.ylabel('$y [kpc]$')
#	plt.savefig('./gif7/f' + "%03d" % (i,) + '.png')
#




#######################################################
#######################################################
#these are not integrating nicely due to the the strong singularities at r=rperi, r=rapo (extra factor of 1/r**2 compared to t_r)
	#
	#def dpsi_func(r, e, l):
	#	#takes peri and apo in m
	#	out = 1./( (r**2.) * np.sqrt( 2.*( e - halo_pot(r/mperkpc) ) - ( l/(r) )**2. ) )
	#	return out
	#
	#def delta_psi(e, l, apo = -1, peri= -1):
	#	#takes e and L in J and J*s, returns the radial period in Gyr
	#	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri(e,l)
	#	part_d_psi, err = scipy.integrate.quad(dpsi_func, peri*(0.0000001)*mperkpc, apo*(0.999999)*mperkpc, args = (e,l))
	#	return 2*l*part_d_psi
	#
	#def azi_period(e, l, apo = -1, peri= -1):
	#	#returns t_psi in Gyr
	#	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri(e,l)
	#	t_psi = 2.*np.pi*radial_period(e,l, apo=apo, peri=peri)/(abs(delta_psi(e,l, apo=apo, peri=peri)))
	#	return t_psi
	#
	#def omega_precess(e, l, apo=-1, peri=-1):
	#	# in rad per Gyr
	#	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri(e,l)
	#	omega_p = (delta_psi(e,l, apo=apo, peri=peri) - 2*np.pi)/radial_period(e,l,apo=apo, peri=peri)
	#	return omega_p
	###########################################################

