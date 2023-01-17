import scvelo as scv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import subprocess

### Settings
cm = 1/2.54
scv.set_figure_params('scvelo')


### malignant cell barcodes
labelFile = "/home/sangho/2020_newbuild/projects/crc_scRNA-seq/4.publication_annotation/malignant_gencode34_v2.1/newvari_v5_gMean/tmp/umapcoord.res0.3.txt"
openLabel = open(labelFile, 'r')
labelLines = openLabel.readlines()
openLabel.close()

header = labelLines[0]
header = header.rstrip().split('\t')
clstIdx = header.index('res0.3'); clstIdx #3
labelLines = labelLines[1:]

clst = [x.split()[clstIdx] for x in labelLines]
umapcoord = [[float(x.split()[1]), float(x.split()[2])] for x in labelLines]

useBcs = [x.split('\t')[0] for x in labelLines]
useBcs = ['_'.join(x.split('-')[2:]) + ':' + x.split('-')[0].replace('-', '_') + 'x' for x in useBcs]
print (len(useBcs)) # 17276

bcClstDict = dict()
for i in range(len(useBcs)):
    bcClstDict[useBcs[i]] = [umapcoord[i], clst[i]]
bcClstDict[list(bcClstDict.keys())[0]]


### merge loom files
dirPath = '/home/sangho/2020_newbuild/projects/crc_scRNA-seq/1.counts_3.1.0_hg19_gencode34/'
dirs = os.listdir(dirPath)
dirs.sort()
dirs = list(filter(lambda x: '_N_' not in x, dirs))

looms = [dirPath + x + '/velocyto/' + x + '.loom' for x in dirs]

cpCmd = 'cp ' + looms[0] + ' merged.loom'
subprocess.getoutput(cpCmd)

import loompy
merged = loompy.connect("merged.loom")
for loom in looms[1:]:
    merged.add_loom(loom)


### malignant cells
crc = scv.read_loom('merged.loom')
malignantcells = list(filter(lambda x: x in useBcs, list(crc.obs_names))); len(malignantcells) # 16531
crc = crc[malignantcells]; crc
crc.obs['res0.3'] = [bcClstDict[x][1] for x in malignantcells]
crc.obs['res0.3'] = crc.obs['res0.3'].astype("category")
crc

plt.figure(figsize = (12*cm, 5*cm))
scv.pl.proportions(crc, groupby = 'res0.3')
plt.savefig('proportions.pdf', bbox_inches="tight")
plt.clf()


### Preprocess
#scv.pp.filter_and_normalize(crc, min_shared_counts = 20, n_top_genes = 2000)
scv.pp.filter_genes(crc, min_shared_counts = 20)
scv.pp.normalize_per_cell(crc)
scv.pp.filter_genes_dispersion(crc, n_top_genes = 2000, subset = False)
scv.pp.log1p(crc)

scv.pp.moments(crc, n_pcs = 50, n_neighbors = 30)
# which runs: scv.pp.pca(crc); scv.pp.neighbors(crc)

#crc.write('malignant_res0.3.h5ad', compression='gzip')

### Velocity tools
scv.tl.recover_dynamics(crc, n_jobs = 12)
scv.tl.velocity(crc, mode = 'dynamical') # [deterministic|stochastic|dynamical]
scv.tl.velocity_graph(crc)


### UMAP
#scv.tl.umap(crc, min_dist = 0.3)
crc.obsm['X_umap'] = np.array([bcClstDict[x][0] for x in malignantcells], dtype = 'float32')
#crc.write('malignant_res0.3.h5ad', compression='gzip')
#crc = scv.read('malignant_res0.3.h5ad')


clst_cols = ["#E41A1C", "#3A85A8", "#629363", "#C4625D", "#FFC81D", "#BF862B", "#EB7AA9"]
### Plotting 
scv.pl.velocity_embedding(crc, basis = 'umap', color= 'res0.3', palette = clst_cols, figsize = (10*cm, 10*cm)) # cellular level
plt.savefig('umap.pdf', bbox_inches="tight"); plt.clf()

scv.pl.velocity_embedding_grid(crc, basis = 'umap', color= 'res0.3', palette = clst_cols, figsize = (8*cm, 8*cm)) # gridlines
plt.savefig('umap.grid.pdf', bbox_inches="tight"); plt.clf()

scv.pl.velocity_embedding_stream(crc, basis = 'umap', color= 'res0.3', palette = clst_cols, legend_loc = "right", figsize = (10*cm, 10*cm)) # streamlines 
plt.savefig('umap.stream.png', bbox_inches="tight"); plt.clf()


### Latent time
scv.tl.latent_time(crc) # dynamical mode
#crc.write('malignant_res0.3.h5ad', compression='gzip')

scv.pl.scatter(crc, color = 'latent_time', color_map = 'gnuplot', figsize = (8*cm, 8*cm), size = 50)
plt.savefig('umap.latenttime.pdf', bbox_inches="tight"); plt.clf()

#top_genes = crc.var['fit_likelihood'].sort_values(ascending=False).index[:300]
#scv.pl.heatmap(crc, var_names = top_genes, sortby = 'latent_time', col_color = 'res0.3', palette = clst_cols, n_convolve = 100)
#plt.savefig('heatmap.latenttime.topvar.png', bbox_inches="tight"); plt.clf()


### Important genes
scv.tl.rank_velocity_genes(crc, groupby = 'res0.3', min_corr = .3)
df = scv.DataFrame(crc.uns['rank_velocity_genes']['names'])
df.head()


### Speed
scv.tl.velocity_confidence(crc)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(crc, c = keys, cmap = 'coolwarm', perc = [5, 95], figsize = (8*cm, 7*cm))
plt.savefig('scatter.velocity_confidence.pdf', bbox_inches="tight"); plt.clf()
#crc.write('malignant_res0.3.h5ad', compression='gzip')


