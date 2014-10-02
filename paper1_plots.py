import cld_fnc as cf

###########################################################
###########################################################
#precession figure

if 1:

	rcircs = np.arange(150)+0.1 #kpc
	vcircs = np.sqrt(cf.gee*cf.m_encl(rcircs)*cf.sminkg/(rcircs*cf.mperkpc))/1000. #km/s
	lcircs = rcircs*cf.mperkpc*vcircs*1000. #J*s
	ecircs = cf.halo_pot(rcircs) + 0.5*(vcircs*1000.)**2

	ls = (np.arange(len(rcircs))+1.)/len(rcircs)*lcircs[-1]

	t_rs = np.zeros((len(rcircs), len(rcircs)))
	t_as = t_rs*0.
	d_ps = t_rs*0.
	apos = t_rs*0.
	peris = t_rs*0.
	alphas = t_rs*0.
	for i in np.arange(len(rcircs)):
		for j in np.arange(len(ls)):
			if ls[j] < lcircs[i]:
				apo, peri = cf.calc_apo_peri_bisect(ecircs[i], ls[j], rcirc=rcircs[i])
				t_rs[j,i] = cf.radial_period(ecircs[i], ls[j], apo=apo, peri=peri)
				t_as[j,i] = cf.azi_period(ecircs[i], ls[j], apo=apo, peri=peri)
				d_ps[j,i] = cf.delta_psi(ecircs[i], ls[j], apo=apo, peri=peri)
				apos[j,i] = apo
				peris[j,i] = peri
				r_t_onehalf = cf.r_t_onehalf(ecircs[i],ls[j], apo=apo, peri=peri)
				alphas[j,i] = cf.delta_psi(ecircs[i], ls[j], apo=apo, peri=r_t_onehalf)
				print j,i

	np.save('t_rs', t_rs)
	np.save('t_as', t_as)
	np.save('d_ps', d_ps)
	np.save('apos', apos)
	np.save('peris', peris)
	np.save('alphas', alphas )

	t_rs   = np.load('t_rs.npy')
	t_as   = np.load('t_as.npy')
	d_ps   = np.load('d_ps.npy')
	apos   = np.load('apos.npy')
	peris  = np.load('peris.npy')
	alphas = np.load('alphas.npy')


	fig = plt.figure()


	clf()


	#cmap = 'rainbow'
	cmap = cm.jet
	cmap = cm.rainbow
	cmap.set_under(color='k')

	plt.subplot(3,3,1)

	plt.plot(rcircs,lcircs*cf.jstokmkpcpersec, label = 'circular orbits')
	plt.xlim([0,150])
	#plt.ylim([min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec])
	plt.text(20,6e26*cf.jstokmkpcpersec, 'inaccesable')
	plt.text(70,4.6e26*cf.jstokmkpcpersec, 'circular orbits',rotation=45)
	plt.text(80,200, 'eccentric orbits')
	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	#plt.legend(['circular orbits'], loc='lower right')


	#plt.contour((t_rs), interpolation = 'nearest' , origin = 'lower',
	#	extent = [min(rcircs),max(rcircs),min(ls),max(ls)],levels = np.arange(50)/60.+.3)



	plt.subplot(3,3,2, aspect='auto')

	plt.imshow(t_rs, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('Gyr')
	plt.title('Radial period [Gyr]')

	plt.subplot(3,3,3)

	#plt.contour(np.log(use_a), interpolation = 'nearest' , origin = 'lower', 
	#	extent = [0.1,1,0.1,1.1], levels = np.arange(20)/20.)

	plt.imshow(t_as, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('Gyr')
	plt.title('Azimuthal period [Gyr]')

	plt.subplot(3,3,4)
	plt.imshow(d_ps, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 3.1, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('rad/orbit')
	plt.title('Azimuthal precession per orbit [rad')

	plt.subplot(3,3,5)
	plt.imshow(apos, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('kpc')
	plt.title('Apogalacticon distance [kpc]')

	plt.subplot(3,3,6)
	plt.imshow(peris, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('kpc')
	plt.title('Perigalacticon distance [kpc]')

	plt.subplot(3,3,7)
	plt.imshow(alphas*180./np.pi, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('degrees')
	plt.title('Alpha [degrees]')

	fig.subplots_adjust(wspace=0.35)
	fig.subplots_adjust(hspace=0.35)


	fig.tight_layout()


###########################################################
###########################################################
#petal truncation - opening angle (alpha) vs velocity ratio (beta)

if 1:

	rcircs = np.arange(150)+0.1 #kpc
	vcircs = np.sqrt(cf.gee*cf.m_encl(rcircs)*cf.sminkg/(rcircs*cf.mperkpc))/1000. #km/s
	lcircs = rcircs*cf.mperkpc*vcircs*1000. #J*s
	ecircs = cf.halo_pot(rcircs) + 0.5*(vcircs*1000.)**2
	ls = (np.arange(len(rcircs))+1.)/len(rcircs)*lcircs[-1]
	t_rs = np.zeros((len(rcircs), len(rcircs)))
	t_as = t_rs*0.
	d_ps = t_rs*0.
	alphas = t_rs*0.
	betas = t_rs*0.
	for i in np.arange(len(rcircs)):
		for j in np.arange(len(ls)):
			if ls[j] < lcircs[i]:
				apo, peri = cf.calc_apo_peri_bisect(ecircs[i], ls[j], rcirc=rcircs[i])
				t_rs[j,i] = cf.radial_period(ecircs[i], ls[j], apo=apo, peri=peri)
				t_as[j,i] = cf.azi_period(ecircs[i], ls[j], apo=apo, peri=peri)
				d_ps[j,i] = cf.delta_psi(ecircs[i], ls[j], apo=apo, peri=peri)
				alphas[j,i] = cf.delta_psi(ecircs[i], ls[j], apo=apo, peri=rcircs[i])
				betas[j,i] = 2.*np.pi*(peri/apo)**2.

	clf()

	plt.plot(ls[ betas[:,15]>0]/lcircs[15],180./np.pi* betas[ betas[:,15]>0,15], c='g', label='beta, 15')
	plt.plot(ls[alphas[:,15]>0]/lcircs[15],180./np.pi*alphas[alphas[:,15]>0,15], c='g', label='alpha, 15', linestyle='--')

	plt.plot(ls[ betas[:,25]>0]/lcircs[25],180./np.pi* betas[ betas[:,25]>0,25], c='b', label='beta, 25')
	plt.plot(ls[alphas[:,25]>0]/lcircs[25],180./np.pi*alphas[alphas[:,25]>0,25], c='b', label='alpha, 25', linestyle='--')

	plt.plot(ls[ betas[:,45]>0]/lcircs[45],180./np.pi* betas[ betas[:,45]>0,45], c='r', label='beta, 45')
	plt.plot(ls[alphas[:,45]>0]/lcircs[45],180./np.pi*alphas[alphas[:,45]>0,45], c='r', label='alpha, 45', linestyle='--')


	plt.xlabel('L/Lcirc')
	plt.ylabel('Truncation angle [deg]')
	plt.legend(['beta, 15','alpha, 15','beta, 25','alpha, 25','beta, 45','alpha, 45'], loc=2)



##########################################################################
##########################################################################
##########################################################################
## calculate everything for a given nfw halo


if 1:
	l_rel = False
	mhalo = 1.77e11
	zinfall = 4.0
	rcircs = np.arange(25)*3.+.25 #kpc
	vcircs = np.sqrt(cf.gee*cf.m_encl_nfw(mhalo,zinfall,rcircs)*cf.sminkg/(rcircs*cf.mperkpc))/1000. #km/s
	lcircs = rcircs*cf.mperkpc*vcircs*1000. #J*s
	ecircs = cf.nfw_pot(mhalo,zinfall,rcircs) + 0.5*(vcircs*1000.)**2

	#es = (np.arange(len(rcircs))+1.)/len(rcircs)*(ecircs[-1]-ecircs[1])+ecircs[1]
	es=ecircs
	t_rs = np.zeros((len(rcircs), len(rcircs)))
	t_as = t_rs*0.
	d_ps = t_rs*0.
	alphas = t_rs*0.
	deisgood = t_rs*0.
	dlisgood = t_rs*0.
	isgood = t_rs*0.
	ddpdl = t_rs*0.
	dtrde = t_rs*0.
	omegaes = t_rs*0.
	omegals = t_rs*0.
	mus = t_rs*0.
	peris= t_rs*0.
	apos = t_rs*0.
	m = 1e8
	e_s = np.zeros((len(rcircs), len(rcircs)))
	l_s = e_s*0.
	e_s_n = e_s*0.
	l_s_n = e_s*0.
	s = e_s*0.
	eps = 0.0000001
	local_l_start = array([[1.-eps,1.,1+eps],[1.-eps,1.,1+eps],[1.-eps,1.,1+eps]])
	local_e_start = array([[1.+eps,1.+eps,1.+eps],[1.,1.,1.],[1.-eps,1.-eps,1.-eps]])
	local_dps = np.zeros((3,3))
	local_trs = np.zeros((3,3))
	nloop = len(rcircs)*len(ls)
	nloopcounter=0.
	for i in np.arange(len(rcircs)):
		if     l_rel: ls = (np.arange(len(rcircs))+1.)/len(rcircs)*lcircs[i]*.99
		if not l_rel: ls = (np.arange(len(rcircs))+1.)/len(rcircs)*lcircs[-1]
		for j in np.arange(len(ls)):
			if ls[j] < lcircs[i]:
				local_l = local_l_start*ls[j]
				local_e = local_e_start*ecircs[i]
				dlstep=(np.gradient(local_l)[1][1,1])
				destep=(np.gradient(local_e)[0][1,1])

				for iter1 in np.arange(3):
					for iter2 in np.arange(3):
						if ((iter1 ==1) | (iter2 == 1)):
							rcirc = cf.get_rcirc_e_nfw(                   mhalo, zinfall, local_e[iter1,iter2], local_l[iter1,iter2])
							apo, peri = cf.calc_apo_peri_bisect_nfw(      mhalo, zinfall, local_e[iter1,iter2], local_l[iter1,iter2], rcirc=rcirc)
							local_dps[iter1,iter2] = cf.delta_psi_nfw(    mhalo, zinfall, local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
							local_trs[iter1,iter2] = cf.radial_period_nfw(mhalo, zinfall, local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)

				ddpdl[j,i] = np.gradient(local_dps,destep,dlstep)[1][1,1]
				dtrde[j,i] = np.gradient(local_trs,destep,dlstep)[0][1,1]

				apo, peri = cf.calc_apo_peri_bisect_nfw(mhalo,zinfall,es[i], ls[j], rcirc=rcircs[i])
				apos[j,i] = apo
				peris[j,i] = peri
				t_rs[j,i] = cf.radial_period_nfw(mhalo,zinfall,es[i], ls[j], apo=apo, peri=peri)
				#t_a not implemented as nfw yet
				#t_as[j,i] = cf.azi_period(es[i], ls[j], apo=apo, peri=peri)
				d_ps[j,i] = cf.delta_psi_nfw(mhalo,zinfall,es[i], ls[j], apo=apo, peri=peri)
				s[j,i], e_s[j,i], l_s[j,i] = cf.scales_nfw(mhalo,zinfall,m,ls[j],peri)
				e_s_n[j,i] = abs(e_s[j,i]/ecircs[i])
				l_s_n[j,i] = l_s[j,i]/ls[j]
				omegaes[j,i], omegals[j,i], alphas[j,i], mus[j,i] = cf.omegaeomegal(mhalo,m/mhalo,zinfall,es[i],ls[j])
			nloopcounter=nloopcounter+1.
			print nloopcounter/nloop



	np.save('t_rs', t_rs)
	np.save('t_as', t_as)
	np.save('d_ps', d_ps)
	np.save('apos', apos)
	np.save('peris', peris)
	np.save('alphas', alphas )

	t_rs   = np.load('t_rs.npy')
	t_as   = np.load('t_as.npy')
	d_ps   = np.load('d_ps.npy')
	apos   = np.load('apos.npy')
	peris  = np.load('peris.npy')
	alphas = np.load('alphas.npy')


	fig = plt.figure()


	clf()


	#cmap = 'rainbow'
	cmap = cm.jet
	cmap = cm.rainbow
	cmap.set_under(color='k')

	plt.subplot(3,3,1)

	plt.plot(rcircs,lcircs*cf.jstokmkpcpersec, label = 'circular orbits')
	plt.xlim([0,150])
	#plt.ylim([min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec])
	plt.text(20,6e26*cf.jstokmkpcpersec, 'inaccesable')
	plt.text(70,4.6e26*cf.jstokmkpcpersec, 'circular orbits',rotation=45)
	plt.text(80,200, 'eccentric orbits')
	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	#plt.legend(['circular orbits'], loc='lower right')


	#plt.contour((t_rs), interpolation = 'nearest' , origin = 'lower',
	#	extent = [min(rcircs),max(rcircs),min(ls),max(ls)],levels = np.arange(50)/60.+.3)



	plt.subplot(3,3,2, aspect='auto')

	plt.imshow(t_rs, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('Gyr')
	plt.title('Radial period [Gyr]')

#	t_as not implemented
	#plt.subplot(3,3,3)
#
	##plt.contour(np.log(use_a), interpolation = 'nearest' , origin = 'lower', 
	##	extent = [0.1,1,0.1,1.1], levels = np.arange(20)/20.)
#
	#plt.imshow(t_as, interpolation = 'nearest', origin='lower',
	#	extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)
#
	#plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	#plt.ylabel('L [km/s kpc]')
	#cbar = plt.colorbar()
	#cbar.set_label('Gyr')
	#plt.title('Azimuthal period [Gyr]')

	plt.subplot(3,3,4)
	plt.imshow(d_ps, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 3.1, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('rad/orbit')
	plt.title('Azimuthal precession per orbit [rad')

	plt.subplot(3,3,5)
	plt.imshow(apos, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('kpc')
	plt.title('Apogalacticon distance [kpc]')

	plt.subplot(3,3,6)
	plt.imshow(peris, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('kpc')
	plt.title('Perigalacticon distance [kpc]')

	plt.subplot(3,3,7)
	plt.imshow(alphas*180./np.pi, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('degrees')
	plt.title('Alpha [degrees]')


	plt.subplot(3,3,8)
	plt.imshow(omegaes, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('rad')
	plt.title('$\Psi_E$')

#plt.subplot(3,3,9)
#plt.imshow(omegals, interpolation = 'nearest', origin='lower',
#	extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

#plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
#plt.ylabel('L [km/s kpc]')
#cbar = plt.colorbar()
#cbar.set_label('rad')
#plt.title('$\Psi_L$')


	plt.subplot(3,3,9)
	plt.imshow(omegals, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]')
	plt.ylabel('L [km/s kpc]')
	cbar = plt.colorbar()
	cbar.set_label('[1]')
	plt.title('$\mu$')





	#fig.subplots_adjust(wspace=0.35)
	#fig.subplots_adjust(hspace=0.35)
	#fig.tight_layout()


#########################################################
#########################################################
# mu evolution
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=72, Om0=0.25)

omegaes = np.zeros(60)
omegals = np.zeros(60)
alphas = np.zeros(60)
mus = np.zeros(60)

mhalo = 1.77e11
zinfall = 4.0
rcirc = 25 #kpc
vcirc = np.sqrt(cf.gee*cf.m_encl_nfw(mhalo,zinfall,rcirc)*cf.sminkg/(rcirc*cf.mperkpc))/1000. #km/s
lcirc = rcirc*cf.mperkpc*vcirc*1000. #J*s
ecirc = cf.nfw_pot(mhalo,zinfall,rcirc) + 0.5*(vcirc*1000.)**2

for i in np.arange(60):
	omegaes[i], omegals[i], alphas[i], mus[i] = cf.omegaeomegal(mhalo,m/mhalo,i/30.+.01,ecirc,lcirc*.05)
	print omegaes[i], omegals[i], alphas[i], mus[i] 

ts = cosmo.age(0).value - cosmo.age(np.arange(60)/30.+.01).value

plt.plot(ts,omegaes, c='r')
plt.plot(ts,omegals, c='b')
plt.plot(ts,alphas, c='g')
plt.plot(ts,mus, c='k')

##########################################################################
##########################################################################
#### cloudiness evolution vs norbs

#equal spacing of e and l
if 1:
	rcircs = np.arange(100)/2+.1 #kpc
	vcircs = np.sqrt(cf.gee*cf.m_encl(rcircs)*cf.sminkg/(rcircs*cf.mperkpc))/1000. #km/s
	lcircs = rcircs*cf.mperkpc*vcircs*1000. #J*s
	ecircs = cf.halo_pot(rcircs) + 0.5*(vcircs*1000.)**2

	es = (np.arange(len(rcircs))+1.)/len(rcircs)*(ecircs[-1]-ecircs[1])+ecircs[1]
	ls = (np.arange(len(rcircs))+1.)/len(rcircs)*lcircs[-1]
	t_rs = np.zeros((len(rcircs), len(rcircs)))
	t_as = t_rs*0.
	d_ps = t_rs*0.
	alphas = t_rs*0.
	betas = t_rs*0.
	deisgood = t_rs*0.
	dlisgood = t_rs*0.
	isgood = t_rs*0.
	m = 1e8
	e_s = np.zeros((len(rcircs), len(rcircs)))
	l_s = e_s*0.
	e_s_n = e_s*0.
	l_s_n = e_s*0.
	s = e_s*0.

	for i in np.arange(len(rcircs)):
		for j in np.arange(len(ls)):
			rcirc = cf.get_rcirc_e(es[i],ls[j])
			vcirc = np.sqrt(cf.gee*cf.m_encl(rcirc)*cf.sminkg/(rcirc*cf.mperkpc))/1000. #km/s
			lcirc = rcirc*cf.mperkpc*vcirc*1000. #J*s
			if ls[j] < lcirc:
				apo, peri = cf.calc_apo_peri_bisect(es[i], ls[j], rcirc=rcirc)
				t_rs[j,i] = cf.radial_period(es[i], ls[j], apo=apo, peri=peri)
				t_as[j,i] = cf.azi_period(es[i], ls[j], apo=apo, peri=peri)
				d_ps[j,i] = cf.delta_psi(es[i], ls[j], apo=apo, peri=peri)
				s[j,i], e_s[j,i], l_s[j,i] = cf.scales(m, ls[j], peri)
				e_s_n[j,i] = abs(e_s[j,i]/ecircs[i])
				l_s_n[j,i] = l_s[j,i]/ls[j]
				alphas[j,i] = cf.delta_psi(es[i], ls[j], apo=apo, peri=rcirc)
				betas[j,i] = 2.*np.pi*(peri/apo)**2.

	dlstep=np.median(np.gradient(ls))
	destep=np.median(np.gradient(es))
	ddpdl = np.gradient(d_ps, dlstep,destep)[0]
	dtrde = np.gradient(t_rs, dlstep,destep)[1]


	deisgood = t_rs*0.
	dlisgood = t_rs*0.
	esizes = e_s*dtrde/t_rs*d_ps
	lsizes = l_s*ddpdl

	for i in np.arange(len(rcircs)):
		for j in np.arange(len(ls)):
			rcirc = cf.get_rcirc_e(es[i],ls[j])
			vcirc = np.sqrt(cf.gee*cf.m_encl(rcirc)*cf.sminkg/(rcirc*cf.mperkpc))/1000. #km/s
			lcirc = rcirc*cf.mperkpc*vcirc*1000. #J*s
			if ls[j] < lcirc:
				if ((i>0) & (i < (len(rcircs)-1))):
					if ((j>0) & (j < (len(rcircs)-1))):
						if ls[j+1] < lcirc:
							dlisgood[j,i]=True
				rcirc = cf.get_rcirc_e(es[i-1],ls[j])
				vcirc = np.sqrt(cf.gee*cf.m_encl(rcirc)*cf.sminkg/(rcirc*cf.mperkpc))/1000. #km/s
				lcirc = rcirc*cf.mperkpc*vcirc*1000. #J*s
				if ((i>0) & (i < (len(rcircs)-1))):
					if ((j>0) & (j < (len(rcircs)-1))):
						if ls[j] < lcirc:
							deisgood[j,i]=True


	if 0:
		plt.clf()

		cmap = cm.hot
		plt.subplot(1,2,1)
		plt.imshow(ddpdl*dlisgood, interpolation = 'nearest', origin='lower',
			extent = (min(es),max(es),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), 
			aspect = 'auto', cmap=cmap, vmin=min(ddpdl[dlisgood==True]), vmax =max(ddpdl[dlisgood==True]) )
		plt.colorbar()

		plt.subplot(1,2,2)

		plt.imshow(dtrde*deisgood, interpolation = 'nearest', origin='lower',
			extent = (min(es),max(es),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), 
			aspect = 'auto', cmap=cmap, vmin=min(dtrde[deisgood==True]), vmax =max(dtrde[deisgood==True]))
		plt.colorbar()

		clf()

		times = np.arange(100)

		esizes = e_s*dtrde/t_rs*d_ps
		lsizes = l_s*ddpdl

		thise,thisl= 45,1
		esize = esizes[thisl,thise]*times
		esize[esize>alphas[thisl,thise]] = alphas[thisl,thise]
		plt.plot(times,esize*180./np.pi)
		plt.plot(times,lsizes[thisl,thise]*times*180./np.pi)

	clf()
	fig=plt.figure()
	times = [0.5,1,2,6]

	for t in np.arange(len(times)):
		plt.subplot(2,4,t+5)
		this_esize = esizes*times[t]
		this_esize[this_esize>alphas] = alphas[this_esize>alphas]
		this_lsize = lsizes*times[t]*180./np.pi
		this_esize = this_esize*180./np.pi
		plt.imshow(this_lsize/this_esize, interpolation='nearest', origin='lower', 
			extent = (min(es),max(es),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), 
			aspect = 'auto', vmax=1, vmin=0)
		#plt.title('Fraction of times spent outside rcirc; %2.2f kpc < rcirc < %2.2f kpc'%(rcircs[0],rcircs[-1]))
		#plt.ylabel('L [km/s kpc]')
		#plt.xlabel('E [J]')
		plt.text(-1.5e11, 8000, 't = %1.1f Gyr'%times[t])
		plt.text(-1.5e11, 7000, 'm = %1.1e M_odot'%(m/1e7))
		plt.colorbar()


	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(hspace=0)
	plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
	plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)
	plt.setp(fig.axes[4].get_yticklabels(), visible=True)
	plt.setp(fig.axes[4].get_xticklabels(), visible=True)

###############################################
#equal spacing of ecirc and lcirc


if 1:
	rcircs = np.arange(150)/3.+.25 #kpc
	vcircs = np.sqrt(cf.gee*cf.m_encl(rcircs)*cf.sminkg/(rcircs*cf.mperkpc))/1000. #km/s
	lcircs = rcircs*cf.mperkpc*vcircs*1000. #J*s
	ecircs = cf.halo_pot(rcircs) + 0.5*(vcircs*1000.)**2

	#es = (np.arange(len(rcircs))+1.)/len(rcircs)*(ecircs[-1]-ecircs[1])+ecircs[1]
	es=ecircs
	t_rs = np.zeros((len(rcircs), len(rcircs)))
	t_as = t_rs*0.
	d_ps = t_rs*0.
	alphas = t_rs*0.
	betas = t_rs*0.
	deisgood = t_rs*0.
	dlisgood = t_rs*0.
	isgood = t_rs*0.
	ddpdl = t_rs*0.
	dtrde = t_rs*0.
	m = 1e8
	e_s = np.zeros((len(rcircs), len(rcircs)))
	l_s = e_s*0.
	e_s_n = e_s*0.
	l_s_n = e_s*0.
	s = e_s*0.
	eps = 0.0000001
	local_l_start = array([[1.-eps,1.,1+eps],[1.-eps,1.,1+eps],[1.-eps,1.,1+eps]])
	local_e_start = array([[1.+eps,1.+eps,1.+eps],[1.,1.,1.],[1.-eps,1.-eps,1.-eps]])
	local_dps = np.zeros((3,3))
	local_trs = np.zeros((3,3))

	for i in np.arange(len(rcircs)):
		ls = (np.arange(len(rcircs))+1.)/len(rcircs)*lcircs[i]*.99
		for j in np.arange(len(ls)):
			local_l = local_l_start*ls[j]
			local_e = local_e_start*ecircs[i]
			dlstep=(np.gradient(local_l)[1][1,1])
			destep=(np.gradient(local_e)[0][1,1])

			for iter1 in np.arange(3):
				for iter2 in np.arange(3):
					if ((iter1 ==1) | (iter2 == 1)):
						rcirc = cf.get_rcirc_e(local_e[iter1,iter2], local_l[iter1,iter2])
						apo, peri = cf.calc_apo_peri_bisect(      local_e[iter1,iter2], local_l[iter1,iter2], rcirc=rcirc)
						local_dps[iter1,iter2] = cf.delta_psi(    local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
						local_trs[iter1,iter2] = cf.radial_period(local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)

			ddpdl[j,i] = np.gradient(local_dps,destep,dlstep)[1][1,1]
			dtrde[j,i] = np.gradient(local_trs,destep,dlstep)[0][1,1]

			apo, peri = cf.calc_apo_peri_bisect(es[i], ls[j], rcirc=rcircs[i])
			t_rs[j,i] = cf.radial_period(es[i], ls[j], apo=apo, peri=peri)
			t_as[j,i] = cf.azi_period(es[i], ls[j], apo=apo, peri=peri)
			d_ps[j,i] = cf.delta_psi(es[i], ls[j], apo=apo, peri=peri)
			s[j,i], e_s[j,i], l_s[j,i] = cf.scales(m, ls[j], peri)
			e_s_n[j,i] = abs(e_s[j,i]/ecircs[i])
			l_s_n[j,i] = l_s[j,i]/ls[j]
			alphas[j,i] = cf.delta_psi(es[i], ls[j], apo=apo, peri=rcircs[i])
			betas[j,i] = 2.*np.pi*(peri/apo)**2.


	esizes = e_s*dtrde/t_rs*d_ps
	lsizes = l_s*ddpdl


	############cloudiness evolution

	clf()
	fig=plt.figure()
	orbs = [1,2,4,8]
	#cmap = cm.rainbow
	#cmap.set_under(color='k')
	#cmap.set_over(color='w')
	for t in np.arange(len(orbs)):
		plt.subplot(2,2,t+1)
		this_esize = esizes*orbs[t]
		this_esize[this_esize>alphas] = alphas[this_esize>alphas]
		this_lsize = lsizes*orbs[t]*180./np.pi
		this_esize = this_esize*180./np.pi
		plt.imshow((this_lsize/this_esize), interpolation='nearest', origin='lower', 
			extent = (min(rcircs),max(rcircs),0,.99), 
			aspect = 'auto', vmax =1, vmin=0)
		#plt.title('Fraction of time spent outside rcirc; %2.2f kpc < rcirc < %2.2f kpc'%(rcircs[0],rcircs[-1]))
		plt.ylabel('$L/L_{circ}$ [1]')
		plt.xlabel('$r_{circ}$ [kpc]')
		plt.xlim([min(rcircs),max(rcircs)])
		plt.ylim([0,.99])
		plt.text(28,0.9, '%1.1f orbits'%orbs[t],color='w')
		plt.text(28,0.82, 'm = %1.1e $M_{\odot}$'%(m),color='w')
		cbar =plt.colorbar()
		cbar.set_label('$\Omega_L/\Omega_E$')
		#plt.scatter([25,25,25],[0.1,0.5, 0.9], s=50, c='k', edgecolor = 'w')


	#########compare cloudiness to sims
	#l=0.1,0.5,0.9 apo snaps
	#cmap = cm.rainbow
	#cmap.set_under(color='k')
	#cmap.set_over(color='w')
	#clf()
	fig = plt.figure()
	d = cf.particles()


	orb1snap = [7,7,7]
	orb2snap = [14,14,14]
	orb4snap = [27,28,28]
	orb8snap = [53,55,56]
	orbsnaps = [orb1snap,orb2snap,orb4snap,orb8snap]
	orbs = [1,2,4,8]
	pickls = ['0.10', '0.50', '0.90']

	clf()


	for t in np.arange(len(orbs)):
		plt.subplot(4,4,(t+1)*4-3)
		this_esize = esizes*orbs[t]
		this_esize[this_esize>alphas] = alphas[this_esize>alphas]
		this_lsize = lsizes*orbs[t]*180./np.pi
		this_esize = this_esize*180./np.pi
		plt.imshow((this_lsize/this_esize), interpolation='nearest', origin='lower', 
			extent = (min(rcircs),max(rcircs),0,.99), 
			aspect = 'auto', vmax =1, vmin=0)
		print (this_lsize/this_esize)[14,75],(this_lsize/this_esize)[14,75]
		plt.xlim([min(rcircs),max(rcircs)])
		plt.ylim([0,.99])
		plt.text(28,0.9, '%1.1f orbits'%orbs[t],color='w')
		plt.text(28,0.82, 'm = %1.1e $M_{\odot}$'%(m),color='w')
		plt.scatter([25,25,25],[0.1,0.5, 0.9], s=50, c='k', edgecolor = 'w')

		orbsnap = orbsnaps[t]


		for j in np.arange(3):
			plt.subplot(4,4,(t+1)*4-j)
			d.init('../L'+pickls[j]+'/',orbsnap[j])
			plt.scatter(d.x,d.y,alpha = 0.06, s=1, edgecolor='none', color = 'k')
			plt.xlim([-48,48])
			plt.ylim([-48,48])
			plt.text(-40,35,'$L/L_{circ} = $'+ pickls[j], color = 'k')

			#plt.subplot(4,4,(t+1)*4-1)
			#d.init('../'+pickls[j],orbsnap[j])
			#plt.scatter(d.x,d.y,alpha = 0.01, s=1, edgecolor='none', color = 'k')

			#plt.subplot(4,4,(t+1)*4)
			#d.init('../'+pickls[j],orbsnap[j])
			#plt.scatter(d.x,d.y,alpha = 0.01, s=1, edgecolor='none', color = 'k')



	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(hspace=0)
	plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
	plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)


	#precession figure - not as useful as absolute L since d/dL direction isn't one of the axes

	if 0: 
		clf()
	
		#cmap = 'rainbow'
		cmap = cm.jet
		cmap = cm.rainbow
		cmap.set_under(color='k')
	
		#plt.subplot(2,2,1)	
		#plt.plot(rcircs,lcircs, label = 'circular orbits')
		#plt.xlim([0,150.1])
		#plt.text(20,6e26, 'inaccesable')
		#plt.text(70,4.6e26, 'circular orbits',rotation=25)
		#plt.text(120,2e26, 'eccentric orbits')
		#plt.xlabel('$R_{circ}$ [kpc]')
		#plt.ylabel('L/1e26 [J s]')
		#plt.legend(['circular orbits'], loc='lower right')
	
	
		#plt.contour((t_rs), interpolation = 'nearest' , origin = 'lower',
		#	extent = [min(rcircs),max(rcircs),min(ls),max(ls)],levels = np.arange(50)/60.+.3)
	
	
	
		plt.subplot(1,3,1)
	
		plt.imshow(t_rs, interpolation = 'nearest', origin='lower',
			extent = (min(rcircs),max(rcircs), 0, 0.99), vmin = 0.001, aspect = 'auto', cmap=cmap)
	
		plt.ylabel('$L/L_{circ}$ [1]')
		plt.xlabel('$r_{circ}$ [kpc]')
		cbar = plt.colorbar()
		cbar.set_label('Gyr')
		plt.title('Radial period [Gyr] as a f(E,L)')
	
		plt.subplot(1,3,2)
	
		#plt.contour(np.log(use_a), interpolation = 'nearest' , origin = 'lower', 
		#	extent = [0.1,1,0.1,1.1], levels = np.arange(20)/20.)
	
		plt.imshow(t_as, interpolation = 'nearest', origin='lower',
			extent = (min(rcircs),max(rcircs), 0, 0.99), vmin = 0.001, aspect = 'auto', cmap=cmap)
	
		plt.ylabel('$L/L_{circ}$ [1]')
		plt.xlabel('$r_{circ}$ [kpc]')
		cbar = plt.colorbar()
		cbar.set_label('Gyr')
		plt.title('Azimuthal period [Gyr] as a f(E,L)')
	
		plt.subplot(1,3,3)
		plt.imshow(d_ps, interpolation = 'nearest', origin='lower',
			extent = (min(rcircs),max(rcircs), 0, 0.99), vmin = 3.1, aspect = 'auto', cmap=cmap)
		plt.ylabel('$L/L_{circ}$ [1]')
		plt.xlabel('$r_{circ}$ [kpc]')
		cbar = plt.colorbar()
		cbar.set_label('rad/orbit')
		plt.title('Azimuthal precession per orbit [rad] as a f(E,L)')


### 4 panel stream/cloud
if 0:
	fig = plt.figure()
	d = cf.particles()
	dire = '/scratch/hendel/clouds/code/nfw/clouds/M2.5e+08/'
	orb4snap = [27,28,28]
	orb8snap = [53,55,56]
	orbsnaps = [orb4snap,orb8snap]
	orbs = [4,8]
	pickls = ['0.50', '0.90']
	rspots = 74
	lspots = [75,135]
	#rcirc of 150 is 74
	#l/lcirc=.9 @ 14
	#l/lcirc=.5 @ 75
	for t in np.arange(len(orbs)):
		orbsnap = orbsnaps[t]
		for j in np.arange(2):
			this_esize = esizes[lspots[j],rspots]*orbs[t]
			if this_esize>alphas[lspots[j],rspots]: this_esize= alphas[lspots[j],rspots]
			this_lsize = lsizes[lspots[j],rspots]*orbs[t]*180./np.pi
			this_esize = this_esize*180./np.pi
			morph = this_lsize/this_esize
			print this_lsize, this_esize, morph
			plt.subplot(2,2,(t+1)*2-j)
			d.init(dire+'L'+pickls[j]+'/',orbsnap[j])
			plt.scatter(d.x,d.y,alpha = 0.06, s=1, edgecolor='none', color = 'k')
			plt.xlim([-48,48])
			plt.ylim([-48,48])
			plt.text(-40,-42,'$\mathrm{L/L_{circ}} = $'+ pickls[j] + ', ' 
				+ '$\mathrm{N_{orb}} = $' + str(orbs[t]) + ', ' + '$\mathrm{\Sigma} = $' + '%1.2f'%morph, color = 'k')
			if (t == 1) & (j == 1):
				plt.xlabel('x [kpc]')
				plt.ylabel('y [kpc]')
			#plt.subplot(4,4,(t+1)*4-1)
			#d.init('../'+pickls[j],orbsnap[j])
			#plt.scatter(d.x,d.y,alpha = 0.01, s=1, edgecolor='none', color = 'k')
			#plt.subplot(4,4,(t+1)*4)
			#d.init('../'+pickls[j],orbsnap[j])
			#plt.scatter(d.x,d.y,alpha = 0.01, s=1, edgecolor='none', color = 'k')
	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(hspace=0)
	plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
	plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)
	plt.savefig('/scratch/hendel/clouds/code/nfw/clouds/proposal_fig1_2.png',bbox_inches='tight')

##########################################################################
##########################################################################
### alpha + beta overlay

if 1:
	rcirc = 25.
	vcirc = np.sqrt(cf.gee*cf.m_encl(rcirc)*cf.sminkg/(rcirc*cf.mperkpc))/1000. #km/s
	lcirc = rcirc*cf.mperkpc*vcirc*1000. #J*s
	ecirc = cf.halo_pot(rcirc) + 0.5*(vcirc*1000.)**2
	ls = (np.arange(10)+1)/10.*lcirc


	x0 = np.zeros(len(ls))
	y0 = x0*0.
	z0 = x0*0.
	vx0 = x0*0.
	vy0 = x0*0.
	vz0 = x0*0.
	alphas = x0*0.
	r_t_onehalfs = x0*0.
	betas = x0*0.
	for i in np.arange(len(ls)):
		apo, peri = cf.calc_apo_peri_bisect(ecirc, ls[i], rcirc=rcirc)
		x0[i]=apo*cf.mperkpc
		vy0[i] = ls[i]/(apo*cf.mperkpc)
		r_t_onehalfs[i] = cf.r_t_onehalf(ecirc,ls[i])
		alphas[i] = cf.delta_psi(ecirc, ls[i], apo=apo, peri=r_t_onehalfs[i])
		betas[i] = 2.*np.pi*(peri/apo)**2.


	dt = cf.secperyr*2e5
	nsteps = 2500.
	xsave = np.zeros((nsteps*2, len(ls)))
	ysave = np.zeros((nsteps*2, len(ls)))
	zsave = np.zeros((nsteps*2, len(ls)))



	#forward integrate orbit from apo

	x,y,z,vx,vy,vz = x0,y0,z0,vx0,vy0,vz0
	ax,ay,az = cf.get_accel_s(x,y,z,vx,vy,vz)
	vx,vy,vz =cf.start_vel(vx,vy,vz,ax,ay,az,dt)

	for i in np.arange(nsteps):
		#plt.scatter(x/cf.mperkpc,y/cf.mperkpc,s=2)
		xsave[i,:] = x
		ysave[i,:] = y
		zsave[i,:] = z
		x,y,z    = cf.steppos(x,y,z,vx,vy,vz,dt)
		ax,ay,az = cf.get_accel_s(x,y,z,vx,vy,vz)
		vx,vy,vz = cf.stepvel(vx,vy,vz,ax,ay,az,dt)

	#go back to apo and then go backwards

	dt = -abs(dt)
	x,y,z,vx,vy,vz = x0,y0,z0,vx0,vy0,vz0
	ax,ay,az = cf.get_accel_s(x,y,z,vx,vy,vz)
	vx,vy,vz =cf.start_vel(vx,vy,vz,ax,ay,az,dt)

	for i in np.arange(nsteps):
		xsave[i+nsteps,:] = x
		ysave[i+nsteps,:] = y
		zsave[i+nsteps,:] = z
		#plt.scatter(x/cf.mperkpc,y/cf.mperkpc,s=2)
		x,y,z    = cf.steppos(x,y,z,vx,vy,vz,dt)
		ax,ay,az = cf.get_accel_s(x,y,z,vx,vy,vz)
		vx,vy,vz = cf.stepvel(vx,vy,vz,ax,ay,az,dt)

	xsave = xsave/cf.mperkpc
	ysave = ysave/cf.mperkpc
	zsave = zsave/cf.mperkpc

	#plt.xlim([-40,40])
	#plt.ylim([-40,40])
	#plt.scatter([0],[0],color='r',s=20)


	fig = plt.figure()
	clf()

	for i in np.arange(9):
		ax = fig.add_subplot(3,3,i+1)
		#plt.scatter(xsave[:,-1],ysave[:,-1] ,s=2, c = 'gray', edgecolor='none')
		import matplotlib.patches
		e1 = patches.Arc((0,0,), r_t_onehalfs[i]*2.,r_t_onehalfs[i]*2., linewidth=2, fill=False, zorder=2, linestyle='dashed')
		ax.add_patch(e1)
		plt.xlim([-40,40])
		plt.ylim([-40,40])
		plt.scatter([0],[0],color='r',s=20)
		plt.scatter(xsave[:,i],ysave[:,i],s=2, c='k', edgecolor='none')
		plt.text(-30,34,'L/Lcirc  = %2.2f'%(ls[i]/lcirc))

		plt.plot([0,40],[0,  40*np.tan(alphas[i]/2.)], color = 'r', linewidth=2)
		plt.plot([0,40],[0, -40*np.tan(alphas[i]/2.)], color = 'r', linewidth=2)
		#plt.plot([0,40],[0,  40*np.tan( betas[i]/2.)], color = 'r', linewidth=2)
		#plt.plot([0,40],[0, -40*np.tan( betas[i]/2.)], color = 'r', linewidth=2)
		if i == 6:
			plt.xlabel('x [kpc]')
			plt.ylabel('y [kpc]')

	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(hspace=0)
	plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
	plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)
	plt.setp(fig.axes[6].get_yticklabels(), visible=True)
	plt.setp(fig.axes[6].get_xticklabels(), visible=True)
	plt.savefig('/scratch/hendel/clouds/code/nfw/clouds/trunc_angle2.png',bbox_inches='tight')


##########################################################################
##########################################################################
#### t(r>rcirc)

if 1:

	rcircs = np.arange(150)/3.+0.1 #kpc
	vcircs = np.sqrt(cf.gee*cf.m_encl(rcircs)*cf.sminkg/(rcircs*cf.mperkpc))/1000. #km/s
	lcircs = rcircs*cf.mperkpc*vcircs*1000. #J*s
	ecircs = cf.halo_pot(rcircs) + 0.5*(vcircs*1000.)**2
	ls = (np.arange(len(rcircs))+1.)/len(rcircs)*lcircs[-1]
	t_rs = np.zeros((len(rcircs), len(rcircs)))
	t_routs = t_rs*0.
	for i in np.arange(len(rcircs)):
		for j in np.arange(len(ls)):
			if ls[j] < lcircs[i]:
				apo, peri = cf.calc_apo_peri_bisect(ecircs[i], ls[j], rcirc=rcircs[i])
				t_rs[j,i] = cf.radial_period(ecircs[i], ls[j], apo=apo, peri=peri)
				t_routs[j,i] = cf.radial_period(ecircs[i], ls[j], apo=apo, peri=rcircs[i])

	clf()
	plt.imshow(t_routs/t_rs, interpolation='nearest', origin='lower', 
		extent = (min(es),max(es),min(ls)*cf.jstokmkpcpersec,max(ls)*cf.jstokmkpcpersec), 
		aspect = 'auto')
	plt.title('Fraction of times spent outside rcirc; %2.2f kpc < rcirc < %2.2f kpc'%(rcircs[0],rcircs[-1]))
	plt.ylabel('L [km/s kpc]')
	plt.xlabel('E [J]')
	plt.colorbar()


##########################################################################
##########################################################################
#### stream vs cloud e and l spatial distribution


#just position

d9 = cf.particles()
d1 = cf.particles()


d1.init('/scratch/hendel/clouds/code/nfw/clouds/M2.5e+08/L0.10/',30)
d9.init('/scratch/hendel/clouds/code/nfw/clouds/M2.5e+08/L0.90/',30)

fig = plt.figure()

if 0: 
	clf()
	plt.subplot(1,2,2,axisbg='black')
	plt.scatter(d1.x,d1.y, alpha = 0.2, s=2, edgecolor='none', color = 'w')
	plt.xlim([-48,48])
	plt.ylim([-48,48])
	plt.text(-40,35,'$L/L_{circ} = 0.1 $', color = 'w', size = 'x-large')
	plt.text(-40,39,'$M = 2.5e8 M_\odot$', color = 'w', size = 'x-large')
	#cbar =plt.colorbar()
	#cbar.set_label('$\Omega_L/\Omega_E$')


	plt.subplot(1,2,1,axisbg='black')
	plt.scatter(d9.x,d9.y, alpha = 0.2, s=2, edgecolor='none', color = 'w')
	plt.xlim([-48,48])
	plt.ylim([-48,48])
	plt.text(-40,35,'$L/L_{circ} = 0.9 $', color = 'w', size = 'x-large')
	plt.text(-40,39,'$M = 2.5e8 M_\odot$', color = 'w', size = 'x-large')
	#cbar =plt.colorbar()
	#cbar.set_label('$\Omega_L/\Omega_E$')



	clf()

	plt.subplot(1,2,2, axisbg='black')
	plt.scatter(d1.x,d1.y, alpha = 0.2, s=3, edgecolor='none', c=(d1.de/d1.e_s), cmap = cm.rainbow)
	plt.xlim([-48,48])
	plt.ylim([-48,48])
	plt.text(-40,35,'$L/L_{circ} = 0.1 $', color = 'w', size = 'x-large')
	plt.text(-40,39,'$M = 2.5e8 M_\odot$', color = 'w', size = 'x-large')
	cbar =plt.colorbar()
	cbar.set_label('$de/e_s$')


	plt.subplot(1,2,1 ,axisbg='black')
	plt.scatter(d9.x,d9.y, alpha = 0.2, s=3, edgecolor='none', c=(d9.de/d9.e_s), cmap = cm.rainbow)
	plt.xlim([-48,48])
	plt.ylim([-48,48])
	plt.text(-40,35,'$L/L_{circ} = 0.9 $', color = 'w', size = 'x-large')
	plt.text(-40,39,'$M = 2.5e8 M_\odot$', color = 'w', size = 'x-large')
	cbar =plt.colorbar()
	cbar.set_label('$de/e_s$')


	clf()

	plt.subplot(1,2,2, axisbg='black')
	plt.scatter(d1.x,d1.y, alpha = 0.2, s=3, edgecolor='none', c=(d1.dl/d1.l_s), cmap = cm.rainbow, vmin = -1.8, vmax = 1.8)
	plt.xlim([-48,48])
	plt.ylim([-48,48])
	plt.text(-40,35,'$L/L_{circ} = 0.1 $', color = 'w', size = 'x-large')
	plt.text(-40,39,'$M = 2.5e8 M_\odot$', color = 'w', size = 'x-large')
	cbar =plt.colorbar()
	cbar.set_label('$dl/l_s$')


	plt.subplot(1,2,1 ,axisbg='black')
	plt.scatter(d9.x,d9.y, alpha = 0.2, s=3, edgecolor='none', c=(d9.dl/d9.l_s), cmap = cm.rainbow)
	plt.xlim([-48,48])
	plt.ylim([-48,48])
	plt.text(-40,35,'$L/L_{circ} = 0.9 $', color = 'w', size = 'x-large')
	plt.text(-40,39,'$M = 2.5e8 M_\odot$', color = 'w', size = 'x-large')
	cbar =plt.colorbar()
	cbar.set_label('$dl/l_s$')

if 1:

	#black background
	#bgcolor ='black'
	#textcolor = 'w'
	#colormap = cm.rainbow

	#white background
	bgcolor ='w'
	textcolor = 'k'
	ssize = 4
	colormap = cm.spectral

	clf()

	plt.subplot(2,2,2, axisbg=bgcolor)
	plt.scatter(d1.x,d1.y, alpha = 0.2, s=ssize, edgecolor='none', c=(d1.de/d1.e_s), cmap = colormap)
	plt.xlim([-48,48])
	plt.ylim([-48,48])
	plt.text(-40,29,'$\mathrm{E}$', color = textcolor, size = 'x-large')
	plt.text(-40,34,'$\mathrm{L/L_{circ} = 0.1 }$', color = textcolor, size = 'x-large')
	plt.text(-40,39,'$\mathrm{M = 6.5x10^8 M_\odot}$', color = textcolor, size = 'x-large')
	#cbar =plt.colorbar()
	#cbar.set_label('$de/e_s$')


	plt.subplot(2,2,1 ,axisbg=bgcolor)
	plt.scatter(d9.x,d9.y, alpha = 0.2, s=ssize, edgecolor='none', c=(d9.de/d9.e_s), cmap = colormap)
	plt.xlim([-48,48])
	plt.ylim([-48,48])
	plt.text(-40,29,'$\mathrm{E}$', color = textcolor, size = 'x-large')
	plt.text(-40,34,'$\mathrm{L/L_{circ} = 0.9 }$', color = textcolor, size = 'x-large')
	plt.text(-40,39,'$\mathrm{M = 6.5x10^8 M_\odot}$', color = textcolor, size = 'x-large')
	#cbar =plt.colorbar()
	#cbar.set_label('$de/e_s$')


	plt.subplot(2,2,4, axisbg=bgcolor)
	plt.scatter(d1.x,d1.y, alpha = 0.2, s=ssize, edgecolor='none', c=(d1.dl/d1.l_s), cmap = colormap, vmin=-2,vmax=2.)
	plt.xlim([-48,48])
	plt.ylim([-48,48])
	plt.text(-40,29,'$\mathrm{L}$', color = textcolor, size = 'x-large')
	plt.text(-40,34,'$\mathrm{L/L_{circ} = 0.1 }$', color = textcolor, size = 'x-large')
	plt.text(-40,39,'$\mathrm{M = 6.5x10^8 M_\odot}$', color = textcolor, size = 'x-large')
	#cbar =plt.colorbar()
	#cbar.set_label('$de/e_s$')


	plt.subplot(2,2,3 ,axisbg=bgcolor)
	plt.scatter(d9.x,d9.y, alpha = 0.2, s=ssize, edgecolor='none', c=(d9.dl/d9.l_s), cmap = colormap)
	plt.xlim([-48,48])
	plt.ylim([-48,48])
	plt.text(-40,29,'$\mathrm{L}$', color = textcolor, size = 'x-large')
	plt.text(-40,34,'$\mathrm{L/L_{circ} = 0.9 }$', color = textcolor, size = 'x-large')
	plt.text(-40,39,'$\mathrm{M = 6.5x10^8 M_\odot}$', color = textcolor, size = 'x-large')
	#cbar =plt.colorbar()
	#cbar.set_label('$de/e_s$')

	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(hspace=0)
	plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
	plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)
	#plt.setp(fig.axes[10].get_yticklabels(), visible=True)
	#plt.setp(fig.axes[17].get_xticklabels(), visible=True)

	plt.savefig('./e_and_l_wbg.png',bbox_inches='tight')



#####################################
#orbit variation ball

if 0: 
	nbods =100000

	dt = cf.secperyr*2e6

	es = -8.3e10*(1+np.random.normal(scale=0.02, size=nbods))
	ls =  0.7e26*(1+np.random.normal(scale=0.08, size =nbods))

	es1 = es
	ls1 = ls

	x = 25.*cf.mperkpc*np.ones(nbods)
	y = 0.*np.ones(nbods)
	z = 0.*np.ones(nbods)

	vy = ls/(x) # in m/s
	vx = -(2.*(es - cf.halo_pot(x/cf.mperkpc) - 0.5*(vy)**2.))**(0.5)
	vz = 0*np.ones(nbods)

	ax,ay,az = cf.get_accel_s(x,y,z,vx,vy,vz)
	vx,vy,vz =cf.start_vel(vx,vy,vz,ax,ay,az,dt)

	x1,y1,z1,vx1,vy1,vz1,ax1,ay1,az1=x,y,z,vx,vy,vz,ax,ay,az


	es = -8.3e10*(1+np.random.normal(scale=0.02, size=nbods))
	ls =  0.7e26*(1+np.random.normal(scale=0.02, size =nbods))

	es2 = es
	ls2 = ls

	x = 25.*cf.mperkpc*np.ones(nbods)
	y = 0.*np.ones(nbods)
	z = 0.*np.ones(nbods)

	vy = ls/(x) # in m/s
	vx = -(2.*(es - cf.halo_pot(x/cf.mperkpc) - 0.5*(vy)**2.))**(0.5)
	vz = 0*np.ones(nbods)

	ax,ay,az = cf.get_accel_s(x,y,z,vx,vy,vz)
	vx,vy,vz =cf.start_vel(vx,vy,vz,ax,ay,az,dt)

	x2,y2,z2,vx2,vy2,vz2,ax2,ay2,az2=x,y,z,vx,vy,vz,ax,ay,az

	nsteps = 3500

	counter=1
	for i in np.arange(nsteps):

		x1,y1,z1    = cf.steppos(x1,y1,z1,vx1,vy1,vz1,dt)
		ax1,ay1,az1 = cf.get_accel_s(x1,y1,z1,vx1,vy1,vz1)
		vx1,vy1,vz1 = cf.stepvel(vx1,vy1,vz1,ax1,ay1,az1,dt)

		x2,y2,z2    = cf.steppos(x2,y2,z2,vx2,vy2,vz2,dt)
		ax2,ay2,az2 = cf.get_accel_s(x2,y2,z2,vx2,vy2,vz2)
		vx2,vy2,vz2 = cf.stepvel(vx2,vy2,vz2,ax2,ay2,az2,dt)


		if 0:
			if mod(i,5) == 0: 
				print i
				clf()
				plt.subplot(1,2,1,axisbg='black')
				plt.scatter(x_center/cf.mperkpc,y_center/cf.mperkpc, s=3, color = 'k', alpha = 0.1, edgecolor='none')
				plt.scatter(x1/cf.mperkpc,y1/cf.mperkpc, s=3, edgecolor='none', color = 'w',alpha = 0.2)

				plt.scatter(0,0, c='r',edgecolor='none')
				plt.title('dL/dE=10')
				plt.xlim([-70,70])
				plt.ylim([-70,70])

				plt.subplot(1,2,2,axisbg='black')
				plt.scatter(x_center/cf.mperkpc,y_center/cf.mperkpc, s=3, color = 'k', alpha = 0.1, edgecolor='none')
				plt.scatter(x2/cf.mperkpc,y2/cf.mperkpc, s=3, edgecolor='none', color = 'w',alpha = 0.2)

				plt.scatter(0,0, c='r',edgecolor='none')
				plt.title('dL/dE=1')
				plt.xlim([-70,70])
				plt.ylim([-70,70])

				fig.subplots_adjust(wspace=0)
				fig.subplots_adjust(hspace=0)
				plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
				plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)

				fig.subplots_adjust(wspace=0)
				fig.subplots_adjust(hspace=0)
				plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
				plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)

				plt.savefig('/scratch/hendel/clouds/code/movie-none/%03d.png'%(counter),bbox_inches='tight')

				#for e
				clf()
				plt.subplot(1,2,1,axisbg='black')
				plt.scatter(x_center/cf.mperkpc,y_center/cf.mperkpc, s=3, color = 'k', alpha = 0.1, edgecolor='none')
				plt.scatter(x1/cf.mperkpc,y1/cf.mperkpc, s=3, edgecolor='none', c=es1/np.mean(es1), alpha = 0.2)

				plt.scatter(0,0, c='r',edgecolor='none')
				plt.title('dL/dE=10')
				plt.xlim([-70,70])
				plt.ylim([-70,70])

				plt.subplot(1,2,2,axisbg='black')
				plt.scatter(x_center/cf.mperkpc,y_center/cf.mperkpc, s=3, color = 'k', alpha = 0.1, edgecolor='none')
				plt.scatter(x2/cf.mperkpc,y2/cf.mperkpc, s=3, edgecolor='none', c=es2/np.mean(es2), alpha = 0.2)

				plt.scatter(0,0, c='r',edgecolor='none')
				plt.title('dL/dE=1')
				plt.xlim([-70,70])
				plt.ylim([-70,70])


				fig.subplots_adjust(wspace=0)
				fig.subplots_adjust(hspace=0)
				plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
				plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)

				plt.savefig('/scratch/hendel/clouds/code/movie-e/%03d.png'%(counter),bbox_inches='tight')


				#for l
				clf()
				plt.subplot(1,2,1,axisbg='black')
				plt.scatter(x_center/cf.mperkpc,y_center/cf.mperkpc, s=3, color = 'k', alpha = 0.1, edgecolor='none')
				plt.scatter(x1/cf.mperkpc,y1/cf.mperkpc, s=3, edgecolor='none', c=ls1/np.mean(ls1),alpha = 0.2)

				plt.scatter(0,0, c='r',edgecolor='none')
				plt.title('dL/dE=10')
				plt.xlim([-70,70])
				plt.ylim([-70,70])

				plt.subplot(1,2,2,axisbg='black')
				plt.scatter(x_center/cf.mperkpc,y_center/cf.mperkpc, s=3, color = 'k', alpha = 0.1, edgecolor='none')
				plt.scatter(x2/cf.mperkpc,y2/cf.mperkpc, s=3, edgecolor='none', c=ls2/np.mean(ls2),alpha = 0.2)

				plt.scatter(0,0, c='r',edgecolor='none')
				plt.title('dL/dE=1')
				plt.xlim([-70,70])
				plt.ylim([-70,70])

				fig.subplots_adjust(wspace=0)
				fig.subplots_adjust(hspace=0)
				plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
				plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)


				plt.savefig('/scratch/hendel/clouds/code/movie-l/%03d.png'%(counter),bbox_inches='tight')
				counter=counter+1


	#black background
	#bgcolor ='black'
	#textcolor = 'w'
	#colormap = cm.rainbow

	#white background
	bgcolor ='w'
	textcolor = 'k'
	ssize = 4
	alpha = 0.2
	colormap = cm.spectral
	limrange = [-70,70]

	clf()

	plt.subplot(2,2,2, axisbg=bgcolor)
	plt.scatter(x1/cf.mperkpc,y1/cf.mperkpc, alpha = alpha, s=ssize, edgecolor='none', c=(es1), cmap = colormap)
	plt.xlim(limrange)
	plt.ylim(limrange)
	plt.text(-60,49,'$\mathrm{E}$', color = textcolor, size = 'x-large')
	plt.text(-60,56,'$\mathrm{\Delta L/ \Delta E = 4 }$', color = textcolor, size = 'x-large')
	#plt.text(-40,29,'$\mathrm{E}$', color = textcolor, size = 'x-large')
	#plt.text(-40,34,'$\mathrm{L/L_{circ} = 0.1 }$', color = textcolor, size = 'x-large')
	#plt.text(-40,39,'$\mathrm{M = 6.5x10^8 M_\odot}$', color = textcolor, size = 'x-large')
	#cbar =plt.colorbar()
	#cbar.set_label('$de/e_s$')


	plt.subplot(2,2,1 ,axisbg=bgcolor)
	plt.scatter(x2/cf.mperkpc,y2/cf.mperkpc, alpha = alpha, s=ssize, edgecolor='none', c=(es2), cmap = colormap)
	plt.xlim(limrange)
	plt.ylim(limrange)
	plt.text(-60,49,'$\mathrm{E}$', color = textcolor, size = 'x-large')
	plt.text(-60,56,'$\mathrm{\Delta L/ \Delta E = 1 }$', color = textcolor, size = 'x-large')
	#plt.text(-40,29,'$\mathrm{E}$', color = textcolor, size = 'x-large')
	#plt.text(-40,34,'$\mathrm{L/L_{circ} = 0.9 }$', color = textcolor, size = 'x-large')
	#plt.text(-40,39,'$\mathrm{M = 6.5x10^8 M_\odot}$', color = textcolor, size = 'x-large')
	#cbar =plt.colorbar()
	#cbar.set_label('$de/e_s$')


	plt.subplot(2,2,4, axisbg=bgcolor)
	plt.scatter(x1/cf.mperkpc,y1/cf.mperkpc, alpha = alpha, s=ssize, edgecolor='none', c=(ls1), cmap = colormap)
	plt.xlim(limrange)
	plt.ylim(limrange)
	plt.text(-60,49,'$\mathrm{L}$', color = textcolor, size = 'x-large')
	plt.text(-60,56,'$\mathrm{\Delta L/ \Delta E = 4 }$', color = textcolor, size = 'x-large')
	#plt.text(-40,29,'$\mathrm{L}$', color = textcolor, size = 'x-large')
	#plt.text(-40,34,'$\mathrm{L/L_{circ} = 0.1 }$', color = textcolor, size = 'x-large')
	#plt.text(-40,39,'$\mathrm{M = 6.5x10^8 M_\odot}$', color = textcolor, size = 'x-large')
	#cbar =plt.colorbar()
	#cbar.set_label('$de/e_s$')


	plt.subplot(2,2,3 ,axisbg=bgcolor)
	plt.scatter(x2/cf.mperkpc,y2/cf.mperkpc, alpha = alpha, s=ssize, edgecolor='none', c=(ls2), cmap = colormap)
	plt.xlim(limrange)
	plt.ylim(limrange)
	plt.text(-60,49,'$\mathrm{L}$', color = textcolor, size = 'x-large')
	plt.text(-60,56,'$\mathrm{\Delta L/ \Delta E = 1 }$', color = textcolor, size = 'x-large')
	#plt.text(-40,29,'$\mathrm{L}$', color = textcolor, size = 'x-large')
	#plt.text(-40,34,'$\mathrm{L/L_{circ} = 0.9 }$', color = textcolor, size = 'x-large')
	#plt.text(-40,39,'$\mathrm{M = 6.5x10^8 M_\odot}$', color = textcolor, size = 'x-large')
	#cbar =plt.colorbar()
	#cbar.set_label('$de/e_s$')

	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(hspace=0)
	plt.setp([a.get_yticklabels() for a in fig.axes[:]], visible=False)
	plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)
	#plt.setp(fig.axes[10].get_yticklabels(), visible=True)
	#plt.setp(fig.axes[17].get_xticklabels(), visible=True)

	plt.savefig('./ball_e_and_l_wbg.png',bbox_inches='tight')


##########################################################################
##########################################################################
####
# figure 2 for proposal - precession & t_r diagrams



if 1:

	rcircs = np.arange(150)+0.1 #kpc
	vcircs = np.sqrt(cf.gee*cf.m_encl(rcircs)*cf.sminkg/(rcircs*cf.mperkpc))/1000. #km/s
	lcircs = rcircs*cf.mperkpc*vcircs*1000. #J*s
	ecircs = cf.halo_pot(rcircs) + 0.5*(vcircs*1000.)**2

	ls = (np.arange(len(rcircs))+1.)/len(rcircs)*lcircs[-1]

	t_rs = np.zeros((len(rcircs), len(rcircs)))
	t_as = t_rs*0.
	d_ps = t_rs*0.

	for i in np.arange(len(rcircs)):
		for j in np.arange(len(ls)):
			if ls[j] < lcircs[i]:
				apo, peri = cf.calc_apo_peri_bisect(ecircs[i], ls[j], rcirc=rcircs[i])
				t_rs[j,i] = cf.radial_period(ecircs[i], ls[j], apo=apo, peri=peri)
				t_as[j,i] = cf.azi_period(ecircs[i], ls[j], apo=apo, peri=peri)
				d_ps[j,i] = cf.delta_psi(ecircs[i], ls[j], apo=apo, peri=peri)

	clf()
	fig=plt.figure()
	#cmap = 'rainbow'
	cmap = cm.jet
	cmap = cm.rainbow
	cmap.set_under(color='k')

	plt.subplot(1,3,2)

	plt.imshow(t_rs, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec/1000.,max(ls)*cf.jstokmkpcpersec/1000.), vmin = 0.001, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]', size ='x-large')
	plt.ylabel('L [km/s Mpc]', size ='x-large')
	cbar = plt.colorbar()
	cbar.set_label('Gyr', size ='x-large')
	plt.title('$\mathrm{T_r}$ [Gyr]', size ='x-large')

	plt.subplot(1,3,3)
	plt.imshow(d_ps*180./np.pi, interpolation = 'nearest', origin='lower',
		extent = (min(rcircs),max(rcircs),min(ls)*cf.jstokmkpcpersec/1000.,max(ls)*cf.jstokmkpcpersec/1000.), vmin = 3.1*180./np.pi, aspect = 'auto', cmap=cmap)

	plt.xlabel('$\mathrm{R_{circ}}$ [kpc]', size ='x-large')
	plt.ylabel('L [km/s Mpc]', size ='x-large')
	cbar = plt.colorbar()
	cbar.set_label('degrees/orbit', size ='x-large')
	plt.title('$\mathrm{\Psi}$ [deg]', size ='x-large')


	rcirc = 25.
	vcirc = np.sqrt(cf.gee*cf.m_encl(rcirc)*cf.sminkg/(rcirc*cf.mperkpc))/1000. #km/s
	lcirc = rcirc*cf.mperkpc*vcirc*1000. #J*s
	ecirc = cf.halo_pot(rcirc) + 0.5*(vcirc*1000.)**2
	ls = (np.arange(10)+1)/10.*lcirc
	ls = [0.3*lcirc]

	x0 = np.zeros(len(ls))
	y0 = x0*0.
	z0 = x0*0.
	vx0 = x0*0.
	vy0 = x0*0.
	vz0 = x0*0.
	alphas = x0*0.
	betas = x0*0.
	for i in np.arange(len(ls)):
		apo, peri = cf.calc_apo_peri_bisect(ecirc, ls[i], rcirc=rcirc)
		x0[i]=apo*cf.mperkpc
		vy0[i] = ls[i]/(apo*cf.mperkpc)
		alphas[i] = cf.delta_psi(ecirc, ls[i], apo=apo, peri=rcirc)
		betas[i] = 2.*np.pi*(peri/apo)**2.


	dt = cf.secperyr*2e5
	nsteps = 5000.
	xsave = np.zeros((nsteps*2, len(ls)))
	ysave = np.zeros((nsteps*2, len(ls)))
	zsave = np.zeros((nsteps*2, len(ls)))



	#forward integrate orbit from apo

	x,y,z,vx,vy,vz = x0,y0,z0,vx0,vy0,vz0
	ax,ay,az = cf.get_accel_s(x,y,z,vx,vy,vz)
	vx,vy,vz =cf.start_vel(vx,vy,vz,ax,ay,az,dt)

	for i in np.arange(nsteps):
		#plt.scatter(x/cf.mperkpc,y/cf.mperkpc,s=2)
		xsave[i,:] = x
		ysave[i,:] = y
		zsave[i,:] = z
		x,y,z    = cf.steppos(x,y,z,vx,vy,vz,dt)
		ax,ay,az = cf.get_accel_s(x,y,z,vx,vy,vz)
		vx,vy,vz = cf.stepvel(vx,vy,vz,ax,ay,az,dt)

	#go back to apo and then go backwards

	dt = -abs(dt)
	x,y,z,vx,vy,vz = x0,y0,z0,vx0,vy0,vz0
	ax,ay,az = cf.get_accel_s(x,y,z,vx,vy,vz)
	vx,vy,vz =cf.start_vel(vx,vy,vz,ax,ay,az,dt)

	for i in np.arange(nsteps):
		xsave[i+nsteps,:] = x
		ysave[i+nsteps,:] = y
		zsave[i+nsteps,:] = z
		#plt.scatter(x/cf.mperkpc,y/cf.mperkpc,s=2)
		x,y,z    = cf.steppos(x,y,z,vx,vy,vz,dt)
		ax,ay,az = cf.get_accel_s(x,y,z,vx,vy,vz)
		vx,vy,vz = cf.stepvel(vx,vy,vz,ax,ay,az,dt)

	xsave = xsave/cf.mperkpc
	ysave = ysave/cf.mperkpc
	zsave = zsave/cf.mperkpc


	import matplotlib.patches
	ax=fig.add_subplot(1,3,1)
	plt.scatter(xsave[:7500],ysave[:7500], c='k', edgecolor='none', s=1)
	plt.xlim([-40,40])
	plt.ylim([-40,40])
	plt.scatter([0],[0],color='r',s=20)
	plt.plot([0,xsave[2743]],[0,ysave[2743]], color = 'b', linewidth=3)
	plt.plot([0,xsave[0]],[0,ysave[0]], color = 'b', linewidth=3)
	plt.plot([0,xsave[1372]],[0,ysave[1372]], color = 'r', linewidth=3)


	xcenter=0
	ycenter=0
	width=20
	height=20
	angle = 0
	theta2 = np.arctan2(ysave[2743],xsave[2743])*180./np.pi

	e1 = patches.Arc((xcenter, ycenter), width, height,
             angle=angle, linewidth=2, fill=False, zorder=2, theta2=theta2, linestyle='dashed')
	ax.add_patch(e1)
	plt.xlabel('x [kpc]', size ='x-large')
	plt.ylabel('y [kpc]', size ='x-large')

	ax.annotate("apogalacticon",
            xy=(xsave[2743], ysave[2743]), xycoords='data', size='large',
            xytext=(-35,-30), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc,rad=.2,angleA=90"),
            )

	ax.annotate("perigalacticon",
            xy=(xsave[1372],ysave[1372]), xycoords='data', size='large',
            xytext=(16,11), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc,rad=.2,angleA=90"),
            )

	#plt.text(-5,12,'$\mathrm{\Delta \Psi}$', size='large')

	ax.annotate('$\mathrm{\Psi}$',
            xy=(-2.5,10), xycoords='data', size='x-large',
            xytext=(-7,18), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc,rad=.2,angleA=90"),
            )

	ax.annotate('',
            xy=(10,0), xycoords='data', size='x-large',
            xytext=(10,.001), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc")
            )



	plt.tight_layout()



##################################################################
#######################################
#es/ls check
#for i in np.arange(9)/10.+0.1:
#	sftp.get('/vega/astro/users/dah2154/nfw/clouds/M2.5e+07/L%1.2f/SCFCEN'%i,'/scratch/hendel/clouds/code/nfw/clouds/M2.5e+07/L%1.2f/SCFCEN'%i)

labels = ['','','','','','','','','']
counter = 0
d = cf.particles()
for i in np.arange(9)/10.+0.1:
	d.init('/scratch/hendel/clouds/code/nfw/clouds/M2.5e+07/L%1.2f/'%i,50)
	print d.e_s
	labels[counter] = 'L=%1.2f'%i
	counter = counter+1
	plt.hist((d.te[d.tub>0]-d.e_c)/d.e_s,histtype='step',bins=np.linspace(-10,10,100), label ='L=%1.2f'%i)

plt.legend(labels)


labels = ['','','','','','','','','']
counter = 0
d = cf.particles()
for i in np.arange(2)/10.+0.1:
	d.init('/scratch/hendel/clouds/code/nfw/clouds/M2.5e+07/L%1.2f/'%i,50)
	print d.l_s
	labels[counter] = 'L=%1.2f'%i
	counter = counter+1
	plt.hist((d.l_mag[d.tub>0]-d.l_mag_c[0])/d.l_s,histtype='step',bins=np.linspace(-10,10,100), label ='L=%1.2f'%i)

plt.legend(labels)

e_s_base,l_s_base = 4825397522.91, 1.64841588413e+25

for i in np.arange(9)/10.+0.1:
	plt.subplot(1,2,1)
	d.init('/scratch/hendel/clouds/code/nfw/clouds/M2.5e+09/L%1.2f/'%i,50)
	plt.scatter((d.l_mag[d.tub>0]-d.l_mag_c[0])/l_s_base,(d.te[d.tub>0]-d.e_c)/e_s_base, edgecolor='none', c='b', s=1, alpha=0.1)
	plt.subplot(1,2,2)
	plt.scatter((d.l_mag[d.tub>0]-d.l_mag_c[0])/d.l_s,(d.te[d.tub>0]-d.e_c)/d.e_s, edgecolor='none', c='b', s=1, alpha=0.1)
	plt.xlim([-10,10])
	plt.ylim([-10,10])

for i in np.arange(9)/10.+0.1:
	plt.subplot(1,2,1)
	d.init('/scratch/hendel/clouds/code/nfw/clouds/M2.5e+07/L%1.2f/'%i,50)
	plt.scatter((d.l_mag[d.tub>0]-d.l_mag_c[0])/l_s_base,(d.te[d.tub>0]-d.e_c)/e_s_base, edgecolor='none', c='r', s=1, alpha=0.1)
	plt.subplot(1,2,2)
	plt.scatter((d.l_mag[d.tub>0]-d.l_mag_c[0])/d.l_s,(d.te[d.tub>0]-d.e_c)/d.e_s, edgecolor='none', c='r', s=1, alpha=0.1)


##################################################################
#######################################
#show phi_l/phi_e



##################################################################
##################################################################
# mu check
clf()
plt.hist(mus,          bins=np.linspace(0,5,100),color='k',cumulative=False, histtype='step')
plt.hist(mus[cls==1.], bins=np.linspace(0,5,100),color='b',cumulative=False, histtype='step')
plt.hist(mus[cls==2.], bins=np.linspace(0,5,100),color='r',cumulative=False, histtype='step')
plt.hist(mus[cls==3.], bins=np.linspace(0,5,100),color='g',cumulative=False, histtype='step')

clf()
plt.subplot(1,1,1,aspect='equal')
plt.scatter(omegals[cls==1.]*180./np.pi,omegae1s[cls==1.]*180./np.pi,c='b',edgecolor='none', label='stream')
plt.scatter(omegals[cls==2.]*180./np.pi,omegae1s[cls==2.]*180./np.pi,c='r',edgecolor='none', label='shell')
plt.scatter(omegals[cls==3.]*180./np.pi,omegae1s[cls==3.]*180./np.pi,c='g',edgecolor='none', label='shell')
plt.plot([0,60],[0*1./.275,60*1./.275],c='k')
plt.xlabel('$\Psi_L$ [degrees]')
plt.ylabel('$\Psi_E^1$ [degrees]')
plt.legend(['$\mu = 0.275$', 'stream', 'shell', 'unclassified'],scatterpoints = 1)
plt.xlim([0,65])
plt.ylim([0,100])


clf()

plt.scatter(omegals[cls==1.],mus[cls==1.],c='b')
plt.scatter(omegals[cls==2.],mus[cls==2.],c='r')
plt.plot([0,2],[0*1./.275,2*1./.275])
plt.xlim([0,2])
plt.ylim([0,2])


dout = np.loadtxt('./classifyer_output5.txt')
clf()
plt.subplot(111,aspect='equal')
plt.scatter(dout[:,5][dout[:,-1]==1]*180./np.pi, dout[:,6][dout[:,-1]==1]*180./np.pi,edgecolor='none', c = 'b')
plt.scatter(dout[:,5][dout[:,-1]==2]*180./np.pi, dout[:,6][dout[:,-1]==2]*180./np.pi,edgecolor='none', c = 'r')
plt.scatter(dout[:,5][dout[:,-1]==3]*180./np.pi, dout[:,6][dout[:,-1]==3]*180./np.pi,edgecolor='none', c = 'g')
plt.scatter(dout[:,5][dout[:,-1]==4]*180./np.pi, dout[:,6][dout[:,-1]==4]*180./np.pi,edgecolor='none', c = 'k')
#plt.plot([0,2],[0*.275,2.*.275])
plt.plot([0,2*180./np.pi],[0,2*180./np.pi])
plt.xlim([0,2*180./np.pi])
plt.ylim([0,2*180./np.pi])
plt.xlabel('$\Psi_E^1$ [degrees]')
plt.ylabel('$\Psi_L$ [degrees]')


dout = np.loadtxt('./classifyer_output5.txt')
clf()
plt.subplot(111,aspect='equal')
plt.scatter(dout[:,-2][dout[:,-1]==1], dout[:,3][dout[:,-1]==1],edgecolor='none', c = 'b')
plt.scatter(dout[:,-2][dout[:,-1]==2], dout[:,3][dout[:,-1]==2],edgecolor='none', c = 'r')
plt.scatter(dout[:,-2][dout[:,-1]==3], dout[:,3][dout[:,-1]==3],edgecolor='none', c = 'g')
plt.scatter(dout[:,-2][dout[:,-1]==4], dout[:,3][dout[:,-1]==4],edgecolor='none', c = 'k')
plt.xlabel('$\Psi_E^1$ [degrees]')
plt.ylabel('$\Psi_L$ [degrees]')
#plt.plot([0,2],[0*.275,2.*.275])
plt.plot([01,01],[0,10])


print 'fraction of streams with mu > 1:', float(np.sum((dout[:,-1]==1) & (dout[:,-2] > 1)))/float(np.sum(dout[:,-1]==1))
print 'fraction of shells  with mu < 1:', float(np.sum((dout[:,-1]==2) & (dout[:,-2] < 1)))/float(np.sum(dout[:,-1]==2))
print 'fraction of debris with mu > 1 not labeled shells ', 1./(float(np.sum((dout[:,-2] > 1)))/float(np.sum((dout[:,-2] > 1) & ((dout[:,-1]==3) | (dout[:,-1]==1) | (dout[:,-1]==3)))))
print 'fraction of debris with mu < 1 not labeled streams', 1./(float(np.sum((dout[:,-2] < 1)))/float(np.sum((dout[:,-2] < 1) & ((dout[:,-1]==3) | (dout[:,-1]==2) | (dout[:,-1]==3)))))




####################################################
#show what type of debris forms in (v_r, v_theta) space


if 0:
	ngrid = 50
	vrs = np.linspace(0.001,1.41, ngrid)
	vts = np.linspace(0.001,1.41, ngrid)
	shout1 = np.zeros((ngrid,ngrid))
	stout1 = np.zeros((ngrid,ngrid))
	shout2 = np.zeros((ngrid,ngrid))
	stout2 = np.zeros((ngrid,ngrid))
	shout3 = np.zeros((ngrid,ngrid))
	stout3 = np.zeros((ngrid,ngrid))
	for i in np.arange(len(vrs)):
		for j in np.arange(len(vts)):
			cf.printclear('%i, %i, %2.2f percent'%(i,j,(100.*(j+i*ngrid))/(ngrid**2.)))
			if vts[j]<np.sqrt(2-vrs[i]**2):
				shout1[i,j], stout1[i,j] = di.integrand_2((0.100,.5,vrs[i], vts[j]), (1e12,'uniform'))
				shout2[i,j], stout2[i,j] = di.integrand_2((0.010,.5,vrs[i], vts[j]), (1e12,'uniform'))
				shout3[i,j], stout3[i,j] = di.integrand_2((0.001,.5,vrs[i], vts[j]), (1e12,'uniform'))

	cf.printclear('done', end=True)

	new2 = shout3/np.max(shout3)*3 +  shout2/np.max(shout2)*2 +  shout1/np.max(shout1)*1 + stout3/np.max(stout)*.01


###
#same but use bisection

if 0:
	from scipy.optimize import bisect
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


	M  = 1e14
	xi = [0.1, 0.01, 0.001, 0.0001] 
	nxi = len(xi)
	z  = 0.5
	vthetamin = 0.005

	from astropy.cosmology import FlatLambdaCDM
	cosmo = FlatLambdaCDM(H0=72, Om0=0.25)


	nelements=50
	vrs = np.linspace(0.001,np.sqrt(2.)-0.001, nelements)
	vthetamax = np.sqrt(2.-vrs**2.)-.0001
	vts1 = np.zeros(nelements)
	vts2 = np.zeros(nelements)
	vts3 = np.zeros(nelements)
	vts4 = np.zeros(nelements)
	vts5 = np.zeros(nelements)
	vtss = [vts1,vts2,vts3,vts4]
	import time
	timeholder=0.
	for i in np.arange(nelements):
		for j in np.arange(nxi):
			if np.sign(ssrc(vthetamin,vrs[i],M,xi[j],z)) != np.sign(ssrc(vthetamax[i],vrs[i],M,xi[j],z)):
				start = time.time()
				vtss[j][i] = bisect(ssrc, vthetamin, vthetamax[i], args = (vrs[i],M,xi[j],z), maxiter =100, disp=False)
				timeholder = timeholder+(time.time()-start)
				cf.printclear('%03d, %03d, %2.2f percent'%(i,j,100.*(i*nxi+j)/(nxi*float(nelements))))

	cf.printclear('%f seconds, %f seconds per bisect'%(timeholder, timeholder/(nelements*nxi)), end=True)



	def ssrc(v_theta, v_r, M, xi, z):
		#shell stream return function
		dnmdxidz = 0.0104*(M/10**12.)**0.133 * xi**-1.995 * np.exp((xi/0.00972)**0.263) * (1+z)**0.0993
		rvir = (M*sminkg/(200.*rho_crit*(4./3.)*np.pi))**(1./3.)/mperkpc # virial radius in kpc
		c = di.halo_concentration(M, z)
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
					rcirc                  = cf.get_rcirc_e_nfw(         M,z,local_e[iter1,iter2], local_l[iter1,iter2]) 
					apo, peri              = cf.calc_apo_peri_bisect_nfw(M,z,local_e[iter1,iter2], local_l[iter1,iter2], rcirc=rcirc)
					if iter1==1: local_dps[iter1,iter2] = cf.delta_psi_nfw(           M,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
					if iter2==1: local_trs[iter1,iter2] = cf.radial_period_nfw(       M,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
					if ((iter1 ==1) & (iter2 == 1)):
						r_t_onehalf        = cf.r_t_onehalf_nfw(         M,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=peri)
						alpha              = cf.delta_psi_nfw(           M,z,local_e[iter1,iter2], local_l[iter1,iter2], apo=apo, peri=r_t_onehalf)

		ddpdl = np.gradient(local_dps,destep,dlstep)[1][1,1]
		dtrde = np.gradient(local_trs,destep,dlstep)[0][1,1]
		s, e_s, l_s = cf.scales_nfw(M, z, msat, l, peri)
		t_r = local_trs[1,1]
		dpsi= local_dps[1,1]
		omega_e = e_s*dtrde/t_r*dpsi*(tinteract/t_r)
		omega_e1 = np.min([alpha,omega_e])
		omega_l = l_s*ddpdl*(tinteract/t_r)
		mu = omega_l/(np.min([alpha,omega_e]))

		return mu-1.

vtss1 = np.load('save_vtss_m1e12_z0.1_100.npy')
vtss2 = np.load('save_vtss_m1e12_z0.5_300.npy')
vtss3 = np.load('save_vtss_m1e13_z0.1_100.npy')
vtss4 = np.load('save_vtss_m1e13_z0.5_100.npy')

colors = ['k','b','r','g','orange']
linestyles = ['-','--','-.',':', '--']
for i in np.arange(5):
	plt.plot(np.linspace(0.001,np.sqrt(2.)-0.001,300), vtss2[i], color = colors[i], linestyle=linestyles[i],lw=2)



colors = ['k','b','r','g','orange']
linestyles = ['-','--','-.',':', '--']
for i in np.arange(5):
	plt.plot(np.linspace(0.001,np.sqrt(2.)-0.001,100), vtss4[i], color = colors[i], linestyle=linestyles[i],lw=2)

plt.plot(np.linspace(0,sqrt(2),1000), np.sqrt(2-np.linspace(0,sqrt(2),1000)**2), c='k', lw=4)


vtss1 = np.load('save_vtss_m1e12_z0.1_100.npy')
vtss2 = np.load('save_vtss_m1e12_z0.5_300.npy')
vtss3 = np.load('save_vtss_m1e10_z0.1_50.npy')
vtss4 = np.load('save_vtss_m1e10_z0.5_50.npy')
vtss5 = np.load('save_vtss_m1e14_z0.1_50.npy')
vtss6 = np.load('save_vtss_m1e14_z0.5_50.npy')

plt.subplot(1,2,1, aspect='equal')
plt.plot(np.linspace(0,sqrt(2),1000), np.sqrt(2-np.linspace(0,sqrt(2),1000)**2), c='k', lw=4)
plt.plot(np.linspace(sqrt(2)-.01,sqrt(2),1000), np.sqrt(2-np.linspace(sqrt(2)-.01,sqrt(2),1000)**2), c='k', lw=4)

for i in np.arange(4):
	plt.plot(np.linspace(0,np.sqrt(2)-.001,len(vtss1[i]))[vtss1[i]>0], vtss1[i][vtss1[i]>0], color = 'r', lw=2)
	plt.plot(np.linspace(0,np.sqrt(2)-.001,len(vtss3[i]))[vtss3[i]>0], vtss3[i][vtss3[i]>0], color = 'g', lw=2)
	plt.plot(np.linspace(0,np.sqrt(2)-.001,len(vtss5[i]))[vtss5[i]>0], vtss5[i][vtss5[i]>0], color = 'b', lw=2)

plt.subplot(1,2,2, aspect='equal')
plt.plot(np.linspace(0,sqrt(2),1000), np.sqrt(2-np.linspace(0,sqrt(2),1000)**2), c='k', lw=4)
plt.plot(np.linspace(sqrt(2)-.01,sqrt(2),1000), np.sqrt(2-np.linspace(sqrt(2)-.01,sqrt(2),1000)**2), c='k', lw=4)

for i in np.arange(4):
	plt.plot(np.linspace(0,np.sqrt(2)-.001,len(vtss2[i]))[vtss2[i]>0], vtss2[i][vtss2[i]>0], color = 'r', lw=2)
	plt.plot(np.linspace(0,np.sqrt(2)-.001,len(vtss4[i]))[vtss4[i]>0], vtss4[i][vtss4[i]>0], color = 'g', lw=2)
	plt.plot(np.linspace(0,np.sqrt(2)-.001,len(vtss6[i]))[vtss6[i]>0], vtss6[i][vtss6[i]>0], color = 'b', lw=2)
