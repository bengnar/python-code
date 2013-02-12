def findclusters(hits_ij):
	'''
	returns all of the clusters of Trues in a boolean matrix
	'''
	blobs = []
	while len(hits_ij)>0:

		ij0 = hits_ij.pop()
		queue_ij = [ij0]

		blob = {}
		# initialize blob
		blob['size'] = 1
		blob['ij'] = [ij0]

		while len(queue_ij)>0:
			center_ij = queue_ij.pop()
			for (di, dj) in [[1, 0], [0, 1], [-1, 0], [0, -1]]:
				test_ij = (center_ij[0]+di, center_ij[1]+dj)
				if hits_ij.count(test_ij)>0:
					blob['size'] += 1
					blob['ij'].append(test_ij)
					queue_ij.append(test_ij)
					hits_ij.remove(test_ij)
		blob['ij'].sort(key = lambda x : x[1]) # sort pixels -- earliest to latest
		blobs.append(blob)
	
	blobs.sort(key = lambda x : x['ij'][0][1])
	return blobs