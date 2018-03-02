#!/usr/bin/env julia

using LaTeXStrings
import PyPlot.plt
using NPZ

img = npzread("img.npy")
mask = npzread("mask.npy")

nx, ny, nv = size(img)

# Estimate how many pixels show emission above some threshold
max_pix = maximum(img)
above = (img .> (1e-3 * max_pix))
println("Total pixels traced ", sum(above))
println("Fraction of pixels brighter than max * 1e-3, ", sum(above)/(nx * ny * nv))

fig, ax = plt[:subplots](nrows=2, ncols=nv, figsize=(1.5 * nv, 1.5))
for k=1:nv
    ax[1,k][:imshow](img[:,:,k], origin="lower", interpolation="none")
    ax[2,k][:imshow](mask[:,:,k], origin="lower", interpolation="none")
end
fig[:subplots_adjust](left=0.03, right=0.97, wspace=0.01)
plt[:savefig]("chmaps.png")
