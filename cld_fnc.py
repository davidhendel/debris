import numpy as np
import scipy.optimize
import scipy.integrate


##############################################################
#general variables

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
#rhalo = a
mhalo=mvir/(np.log(cc+1.)-cc/(cc+1.))


##############################################################
#functions to read snap/com files

def read_snap(dir, snap):
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
	potex = data[:,8] * eu
	tub = data[:,9] * tu

	return m, x, y, z, vx, vy, vz, ep, potex, tub


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
	pars = np.loadtxt(dir + 'SCFPAR', dtype='string')
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

def m_encl(r):
	#returns the halo mass enclosed by r (in kpc) in solar masses
	p=r/rhalo
	m_encl=mhalo*(np.log(p+1.)-p/(p+1.))
	return m_encl

def rho(r):
	#return the density at r
	rho0 = mhalo/(4.*np.pi*rhalo**3.)
	p=r/rhalo
	rho = rho0 / (p * (1+p)**2.)
	return rho

def mass_shell(r):
	#return the mass of a shell at r, for calculating enclosed mass (this matches m_encl)
	return rho(r)*r**2.*4.*np.pi

def scales(m_sat, l_sat, r_scale):
	#returns s = r_j/r, l_s = l/l_c, and e_s = gm/r
	#m_sat in msun, r_scale in kpc
	s = (m_sat/(3.*m_encl(r_scale)))**(1./3.)
	e_s = s * gphys2*m_encl(r_scale)/r_scale
	l_s = s*l_sat
	return s, e_s, l_s

def king_scales(m_sat, l_sat, peri): # e_sat):
	#returns s = r_j/r, l_s = l/l_c, and e_s = gm/r
	#m_sat in msun, r_scale in kpc
	w_peri = l_sat/(peri*mperkpc)**2.
	r_lim = (gee*sminkg*m_sat/(w_peri**2 - dsqrphi_drsqr(peri)))**(1./3.)
	s = r_lim/mperkpc/peri
	e_s = s * peri * mperkpc *dphi_dr(peri)
	l_s = s*l_sat
	return s, e_s, l_s


##############################################################
#orbit properties as f(l,e)
def get_rcirc_l_func(r, ecirc, lcirc):
	return lcirc - r*mperkpc*np.sqrt(gphys2*m_encl(r)/r)

def get_rcirc_l(ecirc,lcirc):
	rcirc = scipy.optimize.newton(get_rcirc_l_func, 30., args = (ecirc,lcirc), maxiter=10000, tol = 1e-9)
	return rcirc

def get_rcirc_e_func(r,ecirc, lcirc):
	return ecirc - (0.5*(np.sqrt(gphys2*m_encl(r)/r))**2. + halo_pot(r))

def get_rcirc_e(ecirc,lcirc):
	rcirc = scipy.optimize.newton(get_rcirc_e_func, 30., args = (ecirc,lcirc), maxiter=10000, tol = 1e-9)
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


def tr_func(r, e, l):
	#takes peri and apo in m
	out = 1./np.sqrt( 2.*( e - halo_pot(r/mperkpc) ) - ( l/(r) )**2. )
	return out

def radial_period(e, l, apo = -1, peri= -1):
	#takes e and L in J and J*s, returns the radial period in Gyr
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri(e,l)
	half_t_r, err = scipy.integrate.quad(tr_func, peri*mperkpc, apo*mperkpc, args = (e,l))
	return (2.*half_t_r)/(secperyr * 10**9.)

###########################################################
#actions

def j_r_func(r, e, l):
	#takes peri and apo in m
	out = np.sqrt( 2.*( e - halo_pot(r/mperkpc) ) - ( l/(r) )**2. )
	return out

def j_r(e, l, apo = -1, peri= -1):
	#takes e and L in J and J*s, returns the radial period in Gyr
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri(e,l)
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


def delta_psi(e, l, apo = -1, peri= -1):
	#takes e and L in J and J*s, returns the radial period in Gyr
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri(e,l)
	part_d_psi, err = scipy.integrate.quad(dpsi_func, 1./(peri*mperkpc), 1./(apo*mperkpc), args = (e,l))
	return 2*l*part_d_psi


def azi_period(e, l, apo = -1, peri= -1):
	#returns t_psi in Gyr
	if (apo < 0. and peri < 0): apo, peri = calc_apo_peri(e,l)
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

def get_accel(x,y,z,vx,vy,vz):
	r = np.sqrt(x**2.+y**2.+z**2.)
	ax = -gee*m_encl(r/mperkpc)*sminkg*x/r**3.#halo_pot(r/mperkpc)*x/r**2
	ay = -gee*m_encl(r/mperkpc)*sminkg*y/r**3.#halo_pot(r/mperkpc)*y/r**2
	az = -gee*m_encl(r/mperkpc)*sminkg*z/r**3.#halo_pot(r/mperkpc)*z/r**2

	return ax,ay,az

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
