import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.ticker import NullFormatter, ScalarFormatter
from matplotlib.collections import LineCollection
import new_cmaps as cmaps

#[:,0] line designation / [:,1] percentage blend / [:,2] wavelength / [:,3] broadened depth
B_hr7512 = np.loadtxt('./HR7512_blend_depths_corr.list',dtype=float)
A_68tau = np.loadtxt('./68Tau_blend_depths_corr.list',dtype=float)
F_procyon = np.loadtxt('./Procyon_blend_depths_corr.list',dtype=float)
G_51peg = np.loadtxt('./51Peg_blend_depths_corr.list',dtype=float)
G_sun = np.loadtxt('./Sun_blend_depths_corr.list',dtype=float)
K_arcturus = np.loadtxt('./Arcturus_blend_depths_corr.list',dtype=float)

plot_list = [G_51peg,G_sun,A_68tau,K_arcturus]#[B_hr7512,A_68tau,F_procyon,G_51peg,G_sun,K_arcturus]
plot_label = ['   51Peg','     Sun','  HR7512','Arcturus']

print len(np.where(G_51peg[:,-1] < 0.98)[0])

'''
def plot_distribs(x,y,plot_label):
		#xy3 = np.vstack([x,y])

		#z3 = gaussian_kde(xy3)(xy3)

		#idx3 = z3.argsort()
		#x,y,z3 = x[idx3], y[idx3], z3[idx3]

		nullfmt = NullFormatter()         # no labels

		# definitions for the axes
		left, width = 0.1, 0.65
		bottom, height = 0.1, 0.65
		bottom_h = left_h = left + width + 0.02

		rect_scatter = [left, bottom, width, height]
		rect_histx = [left, bottom_h, width, 0.2]
		rect_histy = [left_h, bottom, 0.2, height]

		# start with a rectangular Figure
		plt.figure(1, figsize=(8, 8))

		axScatter = plt.axes(rect_scatter)
		axHistx = plt.axes(rect_histx)
		axHisty = plt.axes(rect_histy)

		axScatter.plot([0,100],[0.98,0.98], ls="--", c='royalblue',zorder=2,linewidth=2)
		axScatter.plot([10,10],[0,1], ls="--", c='royalblue',zorder=2,linewidth=2)

		axHisty.plot([0,100000],[0.98,0.98], ls="--", c='royalblue',zorder=2,linewidth=2)
		axHistx.plot([10,10],[0,100000], ls="--", c='royalblue',zorder=2,linewidth=2)

		# no labels
		axHistx.xaxis.set_major_formatter(nullfmt)
		axHisty.yaxis.set_major_formatter(nullfmt)

		# the scatter plot:
		axScatter.scatter(x, y,c='k',s=10,edgecolor='')#,cmap=viridis,zorder=1)

		# now determine nice limits by hand:
		xbinwidth = 1
		ybinwidth = 0.02
		xymaxx = np.max([np.max(0), np.max(100)])
		xymaxy = np.max([np.max(0), np.max(1)])
		xlim = (int(xymaxx/xbinwidth) + 1) * xbinwidth
		ylim = (int(xymaxy/ybinwidth) + 1) * ybinwidth

		axScatter.set_xlim((0,100))
		axScatter.set_ylim((0, 1))

		xbins = np.arange(-xlim, xlim + xbinwidth, xbinwidth)
		ybins = np.arange(-ylim, ylim + ybinwidth, ybinwidth)


		axScatter.set_xlabel('blend percentage (%)')
		axScatter.set_ylabel('Depth (normalised)')

		axHistx.set_ylabel('Count')
		axHisty.set_xlabel('Count')

		axHistx.hist(x, bins=xbins,color='black',zorder=1)
		axHisty.hist(y, bins=ybins, orientation='horizontal',color='black',zorder=1)

		axHistx.set_xlim([0,102])
		axHistx.set_ylim([0,500])

		axHisty.set_xlim([0,1000])
		axHisty.set_ylim([0,1])

		#axHistx.set_yscale('log')
		#axHistx.set_yticks([0.1,1,10, 100, 1000,10000,100000])
		#axHistx.get_yaxis().set_major_formatter(ScalarFormatter())

		#axHisty.set_xscale('log')
		#axHisty.set_xticks([0.1,1,10, 100, 1000,10000,100000])
		#axHisty.get_xaxis().set_major_formatter(ScalarFormatter())

		plt.text(250,1.2,plot_label,fontsize=14)
		plt.show()'''


nullfmt = NullFormatter()         # no labels

# definitions for the axes
left, width = 0.05, 0.315
bottom, height = 0.05, 0.315
bottom_h = left_h = left + width + 0.02

rect_scatter_bl = [left, bottom, width, height]
rect_histx_bl = [left, bottom_h, width, 0.1]
rect_histy_bl = [left_h, bottom, 0.1, height]

rect_scatter_br = [left+0.5, bottom, width, height]
rect_histx_br = [left+0.5, bottom_h, width, 0.1]
rect_histy_br = [left_h+0.5, bottom, 0.1, height]

rect_scatter_tl = [left, bottom, width, height*2]
rect_histx_tl = [left, height*2+left*2, width, 0.2]
rect_histy_tl = [left_h, bottom, 0.1, height*2]

rect_scatter_tr = [left+0.5, bottom, width, height*2]
rect_histx_tr = [left+0.5, height*2+left*2, width, 0.2]
rect_histy_tr = [left_h+0.5, bottom, 0.1, height*2]

# start with a rectangular Figure
plt.figure(1, figsize=(16, 8))

def do_plots(scatter,xhist,yhist,x,y,label):
    axScatter = plt.axes(scatter)
    axHistx = plt.axes(xhist)
    axHisty = plt.axes(yhist)

    axScatter.plot([0,102],[0.02,0.02], ls="-", c='royalblue',zorder=2,linewidth=2)
    axScatter.plot([10,10],[0,1], ls="-", c='royalblue',zorder=2,linewidth=2)

    axHisty.plot([0,100000],[0.02,0.02], ls="-", c='royalblue',zorder=2,linewidth=2)
    axHistx.plot([10,10],[0,100000], ls="-", c='royalblue',zorder=2,linewidth=2)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y,c='k',s=10,edgecolor='')#,cmap=viridis,zorder=1)

    # now determine nice limits by hand:
    xbinwidth = 1
    ybinwidth = 0.02
    xymaxx = np.max([np.max(0), np.max(100)])
    xymaxy = np.max([np.max(0), np.max(1)])
    xlim = (int(xymaxx/xbinwidth) + 1) * xbinwidth
    ylim = (int(xymaxy/ybinwidth) + 1) * ybinwidth

    axScatter.set_xlim((0,102))
    axScatter.set_ylim((0, 1))

    xbins = np.arange(-xlim, xlim + xbinwidth, xbinwidth)
    ybins = np.arange(-ylim, ylim + ybinwidth, ybinwidth)


    axScatter.set_xlabel(r'$\Omega_{core}$ (%)')
    axScatter.set_ylabel(r'd')

    axHistx.set_ylabel('Count')
    axHisty.set_xlabel('Count')

    axHistx.hist(x, bins=xbins,color='black',zorder=1)
    axHisty.hist(y, bins=ybins, orientation='horizontal',color='black',zorder=1)

    axHistx.set_xlim([0,100])
    axHistx.set_ylim([0,500])

    axHisty.set_xlim([0,1000])
    axHisty.set_ylim([0,1])

    plt.text(250,1.2,label,fontsize=18)

    for item in ([axScatter.title, axScatter.xaxis.label, axScatter.yaxis.label] +
             axScatter.get_xticklabels() + axScatter.get_yticklabels()):
        item.set_fontsize(18)
    for item in ([axHisty.title, axHisty.xaxis.label, axHisty.yaxis.label] +
             axHisty.get_xticklabels() + axHisty.get_yticklabels()):
        item.set_fontsize(18)
    for item in ([axHistx.title, axHistx.xaxis.label, axHistx.yaxis.label] +
             axHistx.get_xticklabels() + axHistx.get_yticklabels()):
        item.set_fontsize(18)
    for axisaxis in ['top','bottom','left','right']:
      axScatter.spines[axisaxis].set_linewidth(2)
      axHisty.spines[axisaxis].set_linewidth(2)
      axHistx.spines[axisaxis].set_linewidth(2)

y = 1- plot_list[2][np.where(plot_list[2][:,3] <1)[0],3]
x = plot_list[2][np.where(plot_list[2][:,3] < 1)[0],1]
label = '68 Tau'
do_plots(rect_scatter_tr,rect_histx_tr,rect_histy_tr,x,y,label)

y = 1- plot_list[3][np.where(plot_list[3][:,3] <1)[0],3]
x = plot_list[3][np.where(plot_list[3][:,3] < 1)[0],1]
label = 'Arcturus'
do_plots(rect_scatter_tl,rect_histx_tl,rect_histy_tl,x,y,label)

'''y = 1- plot_list[2][np.where(plot_list[2][:,3] <1)[0],3]
x = plot_list[2][np.where(plot_list[2][:,3] < 1)[0],1]
label = 'HR7512'
do_plots(rect_scatter_br,rect_histx_br,rect_histy_br,x,y,label)

y = 1- plot_list[3][np.where(plot_list[3][:,3] <1)[0],3]
x = plot_list[3][np.where(plot_list[3][:,3] < 1)[0],1]
label = 'Arcturus'
do_plots(rect_scatter_bl,rect_histx_bl,rect_histy_bl,x,y,label)'''
#plt.savefig('./AK_stars.png',bbox_inches='tight')
plt.show()









'''def plot_evolution(x,y):

			# the scatter plot:
			axScatter.plot(x, y,'k-',alpha = 0.1)#,cmap=viridis,zorder=1)

		fig = plt.figure()

		axScatter = fig.add_subplot(111)
		axScatter.set_xlim((0,100))
		axScatter.set_ylim((0, 1))
		axScatter.plot([0,100],[0.98,0.98], ls="--", c='black',zorder=2)
		axScatter.plot([10,10],[0,1], ls="--", c='black',zorder=2)

		for value in range(B_hr7512.shape[0]):
		try:
			x = [K_arcturus[value,1],G_51peg[value,1],F_procyon[value,1],A_68tau[value,1],B_hr7512[value,1]]
			y = [K_arcturus[value,3],G_51peg[value,3],F_procyon[value,3],A_68tau[value,3],B_hr7512[value,3]]

			points = np.array([x, y]).T.reshape(-1, 1, 2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)

			lc = LineCollection(segments, cmap=cmaps.magma,
								norm=plt.Normalize(0, 5),alpha=0.1)
			lc.set_array(np.linspace(0, 5, 5))
			#lc.set_linewidth(3)

			plt.gca().add_collection(lc)
		except IndexError:
			pass

		plt.show()'''
'''
#new_GRADED_file = open('./G_graded_designations.list','w')
full_list = []
tot_list = []
for idx,value in enumerate(plot_list):

	val3 = value[np.where(value[:,3] <1)[0],3]
	val1 = value[np.where(value[:,3] < 1)[0],1]
	plot_distribs(val1,val3,plot_label[idx])

	new_list = []
	for row in range(value.shape[0]):
		if value[row,1] < 10:
			if value[row,3] < 0.98:
				new_list.append(row+1)
				tot_list.append(row+1)
	full_list.append(new_list)
for value in full_list:
	print len(value)

print len(set(tot_list))
'''
'''for line in full_list[0]:
	new_GRADED_file.write(str(line) + '\n')
new_GRADED_file.close()'''
'''
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot([0,100],[0.98,0.98], ls="--", c='black',zorder=2)
ax.plot([10,10],[0,1], ls="--", c='black',zorder=2)

for value in full_list[0]:
	if value not in full_list[1]:
		print G_51peg[value-1,1],G_51peg[value-1,3]
		ax.plot([G_51peg[value-1,1],G_sun[value-1,1]],[G_51peg[value-1,3],G_sun[value-1,3]],'ro-')

for value in full_list[1]:
	if value not in full_list[0]:
		print G_51peg[value-1,1],G_51peg[value-1,3]
		ax.plot([G_51peg[value-1,1],G_sun[value-1,1]],[G_51peg[value-1,3],G_sun[value-1,3]],'bo-')
plt.show()
'''
