# Various DSP helper code for doing some animations
# Champ Darabundit 04/06/21

import numpy as np
import scipy.special as special
import scipy.signal as sig
from scipy.special import binom

def geo_mean(x):
	N = len(x)
	return np.power(np.prod(x), (1./N))

def my_bilinear(b, a, fs=1.0, fo=None, c=None):
	'''
	MY_BILINEAR:
	My version of the bilinear transform. Allows for pre-warping.
	Very much based off scipy signal library function.
	will eventually utilize Sierra method performing filters above Nyquist

	PARAMETERS:
	b	-	Numerator of the analog filter function
	a	-	Denominator of the analog filter function
	fs	-	Sample rate, as ordinary frequency
	fo	-	Pre-warping frequency, as ordinary frequency

	RETURNS
	bz	- 	digital numerator coefficients
	az 	- 	digital denominator coefficients
	'''
	a = np.asarray(a)
	b = np.asarray(b)
	fs = float(fs)
	if c == None:	# We can manually define c
		c = 0.0
		if fo == None:
			c = 2.0 * fs
		else:
			wc = 2.0*np.pi*fo
			c = wc/np.tan(wc/(2*fs))

	D = len(a) - 1
	N = len(b) - 1
	artype = a.dtype
	M = max([N, D])
	Np = M
	Dp = M
	bprime = np.zeros(Np + 1, artype)
	aprime = np.zeros(Dp + 1, artype)
	for j in range(Np + 1):
		val = 0.0
		for i in range(N + 1):
			for k in range(i + 1):
				for l in range(M - i + 1):
					if k + l == j:
						val += (special.comb(i, k) * special.comb(M - i, l) * b[N - i] *
								pow(c, i) * (-1) ** k)
		bprime[j] = val
	for j in range(Dp + 1):
		val = 0.0
		for i in range(D + 1):
			for k in range(i + 1):
				for l in range(M - i + 1):
					if k + l == j:
						val += (special.comb(i, k) * special.comb(M - i, l) * a[D - i] *
								pow(c, i) * (-1) ** k)
		aprime[j] = val
	bprime, aprime = bprime/ aprime[0], aprime/aprime[0]
	bz = bprime.tolist()
	az = aprime.tolist()
	return bz, az

def _lp_arr_hlp(M, g1):
	d = np.zeros(M+1)
	for m in range(M+1):
		d[m] = np.power(g1,(M-m))
	return d

def _bs_arr_hlp(M, N, g1, g2, g3):
	'''
	Helper function to create a general LP->BE transform.
	Separate out top factor into another vector b
	Parameters
	----------
	M - Row size
	N - Column size
	g1 - shift parameter 1
	g2 - shift parameter 2
	g3 - mirror parameter

	Returns
	-------
	A - 2*M+1 x N+1 transformation matrix
	v - an N+1 shift parameter vector
	'''
	A = np.zeros( (2*M+1, N+1) )
	v = np.zeros( (1, N+1) )
	for j in range(N+1):
		for i in range(M+1):
			r_ind = j + 2*i
			if r_ind < 2*M+1:
				A[2*M-r_ind,N-j] = (g2 ** j)*binom(M-j, i)*(g3 ** (M-j-i))
				v[0, N-j] = (g1 ** j)
			else:
				break
	return A, v.flatten()

def mat_lp2bs(b, a, wc, bw, beta=None):
	'''
	Apply bandstop transform using matrix representation
	Parameters
	----------
	b, a - original filter coefficients
	wc - center frequency in rads/s
	bw - bandwidth in rads/s
	beta - original beta is necessary

	Returns
	-------
	bp, ap - bandstop transformed coefficients
	'''
	# Still need to do this to get our order/natural frequency
	z, p, k = sig.tf2zpk(b, a)
	M = len(p)

	if beta is None:
		beta = np.abs(geo_mean(p))

	g1 = beta
	g2 = bw
	g3 = wc ** 2

	A, v = _bs_arr_hlp(M, M, g1, g2, g3)
	# Reshape so these are 2-D
	b = np.pad(b, (M + 1 - len(b), 0), 'constant')
	b = np.multiply(v, b)
	a = np.multiply(v, a)
	a = np.reshape(a, (M + 1, 1))
	b = np.reshape(b, (M + 1, 1))
	bp = np.dot(A,b)
	ap = np.dot(A,a)
	bp, ap = sig.normalize(bp.flatten(),ap.flatten())

	return bp, ap


def bell_filter(N, Q, wc, g=1, g0=1):
	'''
	Bell filter designed using bandstop transfomation
	Parameters
	----------
	N - Half order
	Q - Q factor
	wc - center frequency
	g - gain at peak
	g0 - base gain

	Returns
	-------
	bp, ap - bell filter coefficients
	'''
	b, a = sig.butter(N, 1.0, btype='high', analog=True)
	bw = (1.0*wc)/(1.0*Q)
	bp, ap = mat_lp2bs(b, a, wc, bw, beta=1.0)
	bp = g0*bp
	bp = g*ap + np.pad(bp, (len(ap)-len(bp), 0), 'constant')
	return sig.normalize(bp ,ap)

def stabilize_sys(b, a):
	'''
	Stabilize a system and assume minimum phase
	Parameters
	----------
	b, a - coefficients to stabilize

	Returns
	-------
	b, a - stabilized coefficients
	'''
	z, p, k = sig.tf2zpk(b,a)
	zt = z.copy()
	pt = p.copy()
	zt[zt.real > 0] = -1*zt[zt.real > 0]
	pt[pt.real > 0] = -1*pt[pt.real > 0]
	bp, ap = sig.zpk2tf(zt, pt, k)
	return bp, ap

def _ra(f):
	return (12194*12194*f*f*f*f)/( (f*f + 20.6*20.6)*np.sqrt( (f*f + 107.7*107.7)*(f*f + 737.9*737.9))*(f*f + 12194*12194))

def aweight(f):
	'''
	A weighted frequency weights
	Parameters
	----------
	f - frequency vector to compute weights at

	Returns
	-------
	weights
	'''
	return _ra(f)/_ra(1000)

def dB(x):
	return 20.0*np.log10(np.abs(x))

if __name__ == "__main__":
	pass

