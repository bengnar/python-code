import voc
import glob, os
import itertools
import misc
import Spikes; reload(Spikes)
from scipy.stats import nanmean, nanstd
from matplotlib.pyplot import acorr
basedir = '/Volumes/BOB_SAGET/Vocalization/'

fnames = glob.glob(os.path.join(basedir, 'analysis', 'gui_count_output', '*'))

dtype = np.dtype([('gen', 'S2'), ('cage', 'S8'), ('age', 'i4'), ('chunk', 'i4'), ('ncalls', 'i4'), ('autocorr', '10001f4')])
# 
# db = np.array([], dtype = dtype)
# for fname in fnames:
# 	absol, relat = os.path.split(fname)
# 	cage, age, chunk = relat.split('_')
# 	gen = cage[:2]
# 	age = int(age[1:])
# 	chunk = int(chunk.split('.')[0])
# 	
# 	x = np.loadtxt(fname, ndmin = 2)
# 	ncalls = x.size
# 	
# 	db.resize(db.size + 1)
# 	db[-1] = np.array((gen, cage, age, chunk, ncalls), dtype = dtype)
# 

db = np.array([], dtype = dtype)
lags = np.arange(0, 5000, 10)
empty_autocorr = np.empty(10001) * np.nan
for fname in fnames:
	absol, relat = os.path.split(fname)
	cage, age, chunk = relat.split('_')
	gen = cage[:2]
	age = int(age[1:])
	chunk = int(chunk.split('.')[0])
	
	x = np.loadtxt(fname, ndmin = 2)
	ncalls = x.size
	t2 = np.zeros(65000)
	autocorr = empty_autocorr
	if ncalls > 1:
		t1 = np.int32(x[:, 0]*1000)
		t2[t1] = 1
		tmp = acorr(t2, 'same', maxlags = 5000, normed = False, visible = False)
		autocorr = tmp[1]
	else:
		autocorr = empty_autocorr

	db.resize(db.size + 1)
	db[-1] = np.array((gen, cage, age, chunk, ncalls, autocorr), dtype = dtype)

A = nansum(db['autocorr'], 0)
A_ko = nansum(db[db['gen']=='KO']['autocorr'][:, A.size//2:], 0)
A_wt = nansum(db[db['gen']=='WT']['autocorr'][:, A.size//2:], 0)

A_ko = A_ko[1:] / A_ko[0]
A_wt = A_wt[1:] / A_wt[0]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(Spikes.hamming_smoo(A_wt, windlen = 100), 'b')
ax.plot(Spikes.hamming_smoo(A_ko, windlen = 100), 'g')
ax.set_xlim([0, 4000])
ax.set_xlabel('Lags (ms)')
ax.set_ylabel('P(call)')
ax.set_title('Vocalization onset autocorrelation')
ax.legend(['WT', 'KO'])


		
		

# how many days for each genotype?
cages = np.unique(db['cage'])
dayspercage = dict()
for cage in cages:
	db_ = db[db['cage']==cage]
	dayspercage[cage] = np.unique(db_['age']).size



q = polyfit(db['chunk'], db['ncalls'], 3)
qx = np.linspace(0, 90, 1000)
plot(qx, q[0]*qx**2 + q[1]*qx + q[2])

ax, fig = misc.axis_check(None)

analysis.bar_by_indep('ncalls', 'chunk', db[db['gen']=='WT'], ax = ax, color = 'b')
analysis.bar_by_indep('ncalls', 'chunk', db[db['gen']=='KO'], ax = ax, color = 'r')

ax.set_xlabel('Time (min)')
ax.set_ylabel('No. calls')


gen = np.unique(db['gen'])
age = np.array([(8, 16), (16, 24)])

y = np.empty((age.shape[0], gen.size))
yerr = np.empty((age.shape[0], gen.size))
N = np.empty((age.shape[0], gen.size))
x = np.arange(age.size)
fig = plt.figure()
ax = fig.add_subplot(111)
for i in itertools.product(range(len(gen)), range(len(age))):
	db_ = db[np.vstack((db['gen']==gen[i[0]], db['age']>=age[i[1]][0], db['age']<age[i[1]][1], db['chunk']<=20)).all(0)]
	N[i[0], i[1]] = db_.size
	y[i[0], i[1]] = nanmean(db_['ncalls'])
	yerr[i[0], i[1]] = nanstd(db_['ncalls']) / np.sqrt(db_.size)
	print gen[i[0]], age[i[1]], db_.size, nanmean(db_['ncalls'])





misc.bar2(y, yerr = yerr, groups = [age, gen])

