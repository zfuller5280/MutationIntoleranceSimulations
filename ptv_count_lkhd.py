import numpy as np
import matplotlib

import matplotlib.cm as cm
from scipy.stats import gaussian_kde
from collections import Counter
import random
from scipy.ndimage.filters import gaussian_filter


def get_sample(arr, n_iter=None, sample_size=10,
               fast=True):

    if fast:
        # find the index we last sampled from
        start_idx = (n_iter * sample_size) % 10000
        if start_idx + sample_size >= 10000:
            # shuffle array if we have reached the end and repeat again
            np.random.shuffle(arr)

        return arr[start_idx:start_idx+sample_size]
    else:
        return np.random.choice(arr, sample_size, replace=False)

def kde_scipy(x, obs, bandwidth=0.3, **kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    tol = 1
    try:
        kde = gaussian_kde(x, bw_method=bandwidth, **kwargs)
        return kde.integrate_box(obs-tol,obs+tol)
    except:
        return 0.0

def load_sim(infile):
    sim = np.loadtxt(infile)
    allele_freqs = sim[:, 0]

    return allele_freqs

def make_numpy_array():
    dom_coeffs = np.arange(0.0, 1.01, 0.01)
    sel_coeffs = np.arange(0.01, 1.01, 0.01)
    s = (len(sel_coeffs),len(dom_coeffs),1000000)
    sim_array = np.zeros(s)

    for i, sel in enumerate(sel_coeffs):
        for j, dom in enumerate(dom_coeffs):
            if dom == 0.0:
                sim_infile = (("%.2f_0.0.counts.typical.sim")%(sel))
            else:
                sim_infile = (("%.2f_%.2f.counts.typical.sim")%(sel, dom))
            print sim_infile
            allele_freqs = load_sim(sim_infile)
            sim_array[i,j,:] = allele_freqs

            print sel, dom

    np.save("typical_LoF_counts.npy",sim_array)

def plot_lkhds():
    import matplotlib.pyplot as plt
    sim_array = np.load("typical_LoF_counts.npy")
    sim_array = np.ma.masked_array(sim_array,np.isnan(sim_array))
    mean_array = np.squeeze(sim_array.mean(axis=-1,keepdims=1))

    #Change this value to a different observed count to show probability of that PTV count
    obs_count = 3
    lkhds = np.apply_along_axis(kde_scipy, 2, sim_array, obs_count)
    fig, ax = plt.subplots()
    plt.imshow(lkhds, cmap=cm.seismic)

    cbar = plt.colorbar()
    cbar.set_label('log10(Probability)')
    plt.gca().invert_yaxis()
    plt.ylabel('h')
    plt.xlabel('s')
    ax.set_ylim([0,100])
    ax.set_xlim([0,100])

    CS = plt.contour(gaussian_filter(np.log10(mean_array),sigma=2.5, order=0),levels=[-6.5,-5.5,-4.5,-3.5],colors="black")
    plt.clabel(CS, inline=1, fontsize=10)
    plt.xticks([x for x in range(0,120,25)],[float(x)/100 for x in range(0,120,25)])
    plt.yticks([y for y in range(0,120,25)],[float(y)/100 for y in range(0,120,25)])

    plt.savefig('typical_exac_lof.lkhd.eps', format='eps', dpi=1000)
    plt.show()

def main():
    ##We can uncomment lines here to run the functions
    #make_numpy_array()
    #plot_lkhds()
    
if __name__ == '__main__':
    main()
