# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import sys
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def gettingFacets(filename):
    exe = ["./getFacet", filename]
    result = sp.run(exe, capture_output=True, text=True)
    stderr = result.stderr
    lines = stderr.split("\n")
    segs = []
    skip = False
    for i, line in enumerate(lines[:-1]):
        if line.strip() == '':
            skip = False
        elif not skip:
            temp3 = line.split()
            temp4 = lines[i + 1].split()
            r1, z1 = float(temp3[1]), float(temp3[0])
            r2, z2 = float(temp4[1]), float(temp4[0])
            segs.append(((r1, z1), (r2, z2)))
            segs.append(((-r1, z1), (-r2, z2)))
            skip = True
    return segs

def main():
    nGFS = 27510
    ci = int(sys.argv[1])
    rmin, rmax, zmin, zmax = -5.0, 5.0, 0, 5.0
    lw = 1

    folder = 'Video'
    os.makedirs(folder, exist_ok=True)

    name = f"{ci:04d}_EpsForce.dat"

    if os.path.exists(name):
        print(f"File {name} found! New data will be appended to the file")

    for ti in range(nGFS):
        t = 0.01 * ti
        place = f"intermediate/snapshot-{t:05.4f}"
        if not os.path.exists(place):
            print(f"File {place} not found!")
            continue
        
        exe = f"./getForceNbeta {place} {name}"
        sp.run(exe, shell=True)

        img_name = f"{folder}/{int(t*1000):08d}.png"
        if os.path.exists(img_name):
            print(f"{img_name} Image present!")
            continue

        segs = gettingFacets(place)
        if not segs:
            print(f"Problem in the available file {place}")
            continue

        # Plotting
        fig, ax = plt.subplots(figsize=(19.20, 10.80))
        ax.plot([0, 0], [zmin, zmax], '-.', color='grey', linewidth=lw)
        ax.plot([rmin, rmin], [zmin, zmax], '-', color='black', linewidth=lw)
        ax.plot([rmin, rmax], [zmin, zmin], '-', color='black', linewidth=lw)
        ax.plot([rmin, rmax], [zmax, zmax], '-', color='black', linewidth=lw)
        ax.plot([rmax, rmax], [zmin, zmax], '-', color='black', linewidth=lw)

        line_segments = LineCollection(segs, linewidths=4, colors='green', linestyle='solid')
        ax.add_collection(line_segments)
        ax.set_aspect('equal')
        ax.set_xlim(rmin, rmax)
        ax.set_ylim(zmin, zmax)
        ax.set_title(f'$V_0t/R_0$ = {t:.3f}', fontsize=20)
        ax.axis('off')

        plt.savefig(img_name, bbox_inches="tight")
        plt.close()

        print(f"Done {ti + 1} of {nGFS}")

if __name__ == "__main__":
    main()
