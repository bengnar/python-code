import numpy as np
import os
import sys; sys.path.append('/Users/robert/Documents/pythonwork/A1')
import A1_box, A1_db
import load_session
import load_RF
(params, groups, sessInfos, freqIDmap, basedir, tempdir, sesskeys, groupkeys) = A1_box.load()
DB = A1_db.load_sessions()

fig = plt.figure();
for db in DB:
	fig.clf();
	ax = fig.add_subplot(111);
	ax.imshow(db['zrf'], interpolation = 'nearest', aspect = 'auto');
	ax.axvline(db['bf_man'], color = 'r')
	fig.savefig(os.path.join(tempdir, 'qq_%3.3i' % db['unit']))
		

# present randomized RFs for manual BF manking
fig = plt.figure();
ax = fig.add_subplot(111);
fpath = os.path.join(basedir, 'best_frequency.txt')
rix = np.random.permutation(np.arange(DB.size))
for r, db in enumerate(DB[rix]):
	print r
	fid = open(fpath, 'a')
	zrf_  = db['zrf'].copy()
	zrf_[zrf_<10] = 0
	ax.imshow(zrf_, interpolation = 'nearest', aspect = 'auto');
	plt.show();
	xy = plt.ginput(1, timeout = 0);
	if len(xy) == 0:
		bf = np.nan
		thr = np.nan
	else:
		bf = xy[0][0]
		thr = xy[0][1]
	fid.write('%s\t%3.3i\t%3.3f\t%3.3f\n' % (db['sess'], db['unit'], bf, thr))
	fid.close()


# apply manual BFs to DB
bf_txt = np.loadtxt(fpath, 'str')
for db in DB:
	ix = (bf_txt[:, 0] == db['sess']) & (np.int32(bf_txt[:, 1]) == db['unit'])
	db['bf_man'] = np.float32(bf_txt[ix, 2])
	db['thr_man'] = np.float32(bf_txt[ix, 3])
	print '%3.3f\t%s' % (db['bf_man'], bf_txt[ix, 2][0])

for db in DB:
	if (db['sess'] == 'exp20110727') & (db['unit'] == 68):
		db['goodunit'] = True
	
for db in DB
fig = plt.figure();
ax = fig.add_subplot(111)
for db in DB:
	ax.cla();
	ax.imshow(db['rf'], interpolation = 'nearest', aspect = 'auto')
	ax.axvline(db['bf_com'], color = 'r')
	ax.axvline(db['bf'], color = 'y')
	fig.savefig(os.path.join(basedir, db['sess'], 'analysis/sheets/', 'cc_%3.3i.png' % db['unit']))

DB_ = DB[DB['sess'] == sesskey]

fig = plt.figure()
ax = fig.add_subplot(111)
for db in DB_:
	ax.cla();
	rf = db['zrf']
	ax.imshow(db['rf'], interpolation = 'nearest', aspect = 'auto')
	ax.axvline(db['bf_com'], color = 'r')
	fig.savefig(os.path.join(tempdir, 'cc_%3.3i.png' % db['unit']))

# save out the maps
analyze = 'bf_man'
fig = plt.figure();
# DB2 = np.empty(0, dtype = DB.dtype.descr+[('area', 'f4')])
for sesskey in sesskeys:
	fig.clf();
	ax = fig.add_subplot(111);
	DB_ = DB[DB['sess'] == sesskey] # make a sub-DB with only this session's penetrations
	try:
		DB_ = A1_db.add_field(DB_, [('area', 'f4')])
	except:
		pass

	border = np.loadtxt(os.path.join(basedir, sesskey, 'border.txt'), 'float32') # load the border coordinates
	cf = DB_[analyze] # variable of CFs
	cf[DB_['goodunit']==False] = np.nan # set the bad units cf to nan
	coords = np.concatenate((DB_['coord'][:, 0, :], border), 0) # the penetration AND border coordinates
	cfs = np.hstack((cf, np.ndarray(border.shape[0])*np.nan)) # the pentration AND border CFs
	x, y, a = load_session.load_session.calc_voronoi(coords)
	load_session.load_session.plot_voronoi(cfs, coords, x, y, unitIDs = DB_['unit'], alpha = 0.75, outline_range = np.array([27, 33]))
	ax.scatter(border[:, 0], border[:, 1], marker = 'x')
	ax.set_xlim((border[:, 0].min(), border[:, 0].max()))
	ax.set_ylim((border[:, 1].min(), border[:, 1].max()))
	ax.set_title('%s -- %s' % (analyze, sesskey))
	plt.show();
	fig.savefig(os.path.join(basedir, sesskey, 'analysis', 'map_%s_%s.png' % (analyze, sesskey)))
	for i, db in enumerate(DB_):
		DB_[i]['area'] = a[i]

	A1_db.save(DB_)
	# DB2 = np.append(DB2, DB_)


# calculate percent A1 tuned to different frequencies
fig = plt.figure();
ax = fig.add_subplot(111);
analyze = 'bf_man'
fstep = 6
fbins = np.arange(0, 55, fstep)
clrs = {'exp' : 'r', 'ctl' : 'b'}
for i, sesskey in enumerate(sesskeys):
	print sesskey
	DB_ = DB[(DB['sess'] == sesskey) & (DB['goodunit'])]
	cf = DB_[analyze]
	area = DB_['area']
	n = np.array([])
	for fbin in fbins[:-1]:
		n = np.append(n, np.sum(area[((fbin-(fstep/2))<cf) & (cf<(fbin+(fstep/2)))]))
	n = n/np.sum(area)
	ax.bar(fbins[:-1]+0.8*i, n, color = clrs[DB_['group'][0]]);

ax.legend(sesskeys);
plt.show();
fig.savefig('a1_bf_dist_%s.png' % analyze)

fig = plt.figure();
ax = fig.add_subplot(111);
ix = DB['sess'] == sesskey
ax.scatter(DB[ix]['coord'][:, 0, 0], DB[ix]['coord'][:, 0, 1])
border = np.loadtxt(os.path.join(SESS.basedir, 'border.txt'), 'float32')
coords = np.concatenate((DB['coord'][:, 0, :], border), 0)
cfs = np.hstack((DB['bf_man'], np.ndarray(border.shape[0])*np.nan))

x, y, a = load_session.load_session.calc_voronoi(coords)
load_session.load_session.plot_voronoi(cfs, coords, x, y)
plt.scatter(border[:, 0], border[:, 1], marker = 'x')

x = np.array(x)
y = np.array(y)

a = y[z[0], 0]
b = y[z[0], 1]
c = y[z[0], 2]

for Z in z:
	plot(x[Z[1]], y[Z[2]])
	
for i, (r, g, b) in enumerate(zip(nr, ng, nb)):
	scatter(i, i, color = [r, g, b])
	
x
omit = np.array([65, 39, 76, 60, 91, 93, 78, 81])

sx = set(x)

fig = plt.figure();
ax1 = fig.add_subplot(211);
ax2 = fig.add_subplot(212);
for rf in DB['rf']:
	ax1.cla(); ax2.cla();
	ax1.imshow(rf)
	data = np.sum(rf, 0)
	ax2.plot(np.sum(rf, 0))
	

gaussian = lambda x: 3*exp(-(30-x)**2/20.)

data = gaussian(np.arange(100))

plt.plot(data)

X = np.arange(data.size)
x = np.sum(X*data)/np.sum(data)
width = np.sqrt(np.abs(np.sum((X-x)**2*data)/np.sum(data)))

fit = lambda t : data.max()*exp(-(t-x)**2/(2*width**2))

plot(fit(X))

show()

RF[RF < ]
weight = np.tile(np.linspace(10, 13, 8), (51, 1)).T
RF_weighted = RF * weight
imshow(RF_weighted)

fig = plt.figure();
ax1 = fig.add_subplot(221);
ax2 = fig.add_subplot(222);
ax3 = fig.add_subplot(223);
ax4 = fig.add_subplot(224);
for db in DB:
	print db['unit']
	ax1.cla(); ax2.cla(); ax3.cla(); ax4.cla();
	ax1.cla(); ax2.cla();
	rf = db['rf']
	rf2 = rf.copy()
	rf2[rf<0.2*np.max(rf)] = 0
	rf3 = rf.copy()
	rf3 = (rf3 - db['baselinefr']) / db['baselinestd']
	rf3[rf3<25] = 0
	rf4 = rf3.copy()
	try:
		rf4 = load_RF.load_RF.findmaxcluster(rf4)
	except:
		rf4 = np.nan((2, 2))
	ax1.imshow(rf, aspect = 'auto', interpolation = 'nearest')
	ax2.imshow(rf2, aspect = 'auto', interpolation = 'nearest')
	ax3.imshow(rf3, aspect = 'auto', interpolation = 'nearest')
	ax4.imshow(rf4, aspect = 'auto', interpolation = 'nearest')
	plt.show();
	fig.suptitle(np.str(db['unit']))
	fig.savefig(os.path.join(tempdir, 'aa_%3.3i.png' % db['unit']))



