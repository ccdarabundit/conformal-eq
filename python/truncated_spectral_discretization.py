# TRUNCATED SPECTRAL DISCRETIZATION - CHAMP DARABUNDIT
import numpy as np
import scipy.optimize
from helpers import *
from helpers import _bs_arr_hlp, _lp_arr_hlp


def optimize_g(fs, N = 2**12, weighted=False):
	'''
	Optimize the g-parameter
	Parameters
	----------
	fs - sampling rate in Hz
	N - optional, number of points optimize over
	weighted - using weighting or not

	Returns
	-------

	'''
	w = np.array(range(1, N + 1)) / N * (np.pi * fs)
	f = w / (2.0 * np.pi)

	def peakmap(w):
		wo = np.pi * fs
		return (2 * wo * wo * w) / (wo * wo + w * w)

	def inv_bs(w, g):
		wo = np.pi*fs
		return -(g - np.sqrt(g*g + 4*w*w*wo*wo))/(2*w)

	def bilinear_map(w):
		c = 2.0*fs
		return c*np.tan(w/c)

	def peakmap_full(w, g):
		return peakmap(inv_bs(bilinear_map(w), g))

	g = 2.0*np.pi*fs*np.sqrt((np.pi*fs)**2 - (np.pi*fs/2)**2)

	if weighted:
		weight = aweight(f)

		def peakmap_full_weight(w, g):
			map = peakmap_full(w, g)
			return np.multiply(weight, map)

		return scipy.optimize.curve_fit(peakmap_full_weight, w, np.multiply(weight, w), p0=[g])[0]
	else:
		return scipy.optimize.curve_fit(peakmap_full, w, w, p0=[g])[0]


def tsd(b, a, fs, beta=None, wo=None, g=None, dig=True):
	'''
	Truncated spectral discretization method
	Parameters
	----------
	b - numerator coefficients
	a - denominator coefficients
	fs - sampling rate (Hz)
	beta - optional, warping frequency
	wo - optional, mirror frequency
	g - optional, g-parameter

	Returns
	-------
	bp, ap - corrected numerator and denominator coeffcients
	'''
	p = np.roots(a)

	if wo is None:
		wo = np.pi*fs 	# This is the mirror point, fs in radians

	if beta is None:	# Use pole locations to determine beta
		beta = np.real(geo_mean(p))

	if beta > wo:
		# Mirror beta
		beta = (wo*wo)/beta

	M = len(a) - 1
	F = _lp_arr_hlp(M, wo)

	# 1. Forward transform
	# g1 = -2*wc, g2 = wo, g3 = wo**2
	Dp, G2 = _bs_arr_hlp(M, M, g1=2.0, g2=-1.0, g3=-1.0)
	b = np.pad(b, (M + 1 - len(b), 0), 'constant')
	# Scale by F
	bp = Dp@np.diag(G2)@np.diag(F)@b.reshape(M+1, 1)
	ap = Dp@np.diag(G2)@np.diag(F)@a.reshape(M+1, 1)

	# 2. Stabilize
	bp, ap = stabilize_sys(bp.flatten(), ap.flatten())

	# Compute g parameter
	bnorm = beta / wo
	c = 2*fs
	bh = (c * np.tan(beta / c))/wo
	# Calculate g
	if g == None:
		g = (2*bh*(1 - bnorm**2)**(1/2))/bnorm
	else:
		g = g/(wo*wo)

	# Compute inverse transform
	Ds, G = _bs_arr_hlp(M, M, g1=g, g2=1.0, g3=1.0)
	DsInv = np.linalg.pinv(Ds)
	MP = len(ap)

	if len(bp) < MP:
		bp = np.concatenate((np.zeros(MP - len(bp)), bp))

	bp = np.reshape(bp, (MP, 1))
	ap = np.reshape(ap, (MP, 1))

	# 3. Inverse bandstop transform
	bp = np.diag(F**(-1))@np.diag(G**(-1))@DsInv@bp
	ap = np.diag(F**(-1))@np.diag(G**(-1))@DsInv@ap

	bp, ap = sig.normalize(bp.flatten(), ap.flatten())
	if dig:
		return my_bilinear(bp, ap, fs)
	else:
		return bp, ap



if __name__=='__main__':
	import matplotlib.pyplot as plt
	fig, ax = plt.subplots()
	fs = 48e3
	N = 2**13
	f = np.array(range(N))/N*fs
	w = 2.0*np.pi*f
	wz = w/fs
	# Demo with bell filter
	b, a = bell_filter(2, .75, 2.0*np.pi*25e3)
	_, ha = sig.freqs(b, a, worN=w)
	bpz, apz = tsd(b, a, fs)
	_, hz = sig.freqz(bpz, apz, worN=wz)
	# Now do optimization
	g = optimize_g(fs, weighted=True)
	bpz2, apz2 = tsd(b, a, fs, g=g)
	_, hz2 = sig.freqz(bpz2, apz2, worN=wz)
	ax.semilogx(f, dB(ha), label='Continuous-Time')
	ax.semilogx(f, dB(hz), label='TSD')
	ax.semilogx(f, dB(hz2), label='Optimized g')
	ax.set(
		xlabel='Frequency (Hz)',
		ylabel='Magnitude (dB)',
		title='Truncated Spectral Discretization'
	)
	ax.legend(loc='best')
	ax.grid(True)
	plt.show()
