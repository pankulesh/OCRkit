#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import sys
size = sys.argv[1]
fig, axes = plt.subplots(6, 5, figsize=(25, 30), subplot_kw={'xticks': [], 'yticks': []})
fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
img = []
for i in [11, 12, 13, 14, 15, 16]:
    for j in [1, 2, 3, 4, 5]:
        img.append(mpimg.imread('map{0}_{1}_{2}.png'.format(size, i, j)))
for ax, image in zip(axes.flat, img):
    ax.axis('off')
    ax.imshow(image)
fig.savefig('maps{}.png'.format(size), dpi=300, frameon=False)
