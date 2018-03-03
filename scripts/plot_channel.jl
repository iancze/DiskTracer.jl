#!/usr/bin/env julia

using LaTeXStrings
import PyPlot.plt
using NPZ

img = npzread("img.npy")
mask = npzread("mask.npy")
mask_height = npzread("mask_height.npy")

nx, ny, nv = size(img)

# Estimate how many pixels show emission above some threshold
max_pix = maximum(img)
println("Max pix value ", max_pix)
above = (img .> (1e-5 * max_pix))
println("Total pixels traced ", sum(above))
println("Fraction of pixes unmasked ", sum(mask)/(nx * ny * nv))
println("Fraction of pixes unmasked-height ", sum(mask_height)/(nx * ny * nv))
println("Fraction of pixels brighter than max * 1e-3, ", sum(above)/(nx * ny * nv))

println("Total pixels mask ", sum(mask))
println("Total pixels mask height ", sum(mask_height))

fig, ax = plt[:subplots](nrows=3, ncols=nv, figsize=(1.5 * nv, 4.5))
for k=1:nv
    ax[1,k][:imshow](img[:,:,k], origin="lower", interpolation="none")
    ax[2,k][:imshow](mask[:,:,k], origin="lower", interpolation="none")
    ax[3,k][:imshow](mask_height[:,:,k], origin="lower", interpolation="none")
end
fig[:subplots_adjust](left=0.03, right=0.97, wspace=0.01)
plt[:savefig]("chmaps.png")
