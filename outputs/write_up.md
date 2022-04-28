Cleaning
================

This document is a brief summary of the surface cleaning experiment
performed by Beth in the public health labs at the University of Leeds.

# Experimental set up

5 coupons of PVC (3 x 3cm) were inoculated with with 1ml of *S. aureus*
solution at
![1 \\times 10^9](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%20%5Ctimes%2010%5E9 "1 \times 10^9")
cfu/ml and allowed to dry for 60 min in a stable condition (inside a
fume cupboard - would need to measure typical humidity or look to see if
those data exist).

The coupons were then cleaned with one of the following:

-   Disinfectant wipe (Tuffie - 70% solution of Isopropyl Alcohol BP),
-   detergent wipe (Tuffie - Contains amongst other ingredients less
    than 5% cationic surfactants, amphoteric surfactants and EDTA
    (Vernacare, 2018). They contain a mixture of non-ionic constituents
    at neutral pH),
-   distilled water wipe (Tuffie - alcohol wipes, after being dried,
    washed and re-soaked in distilled water.)
-   no cleaning - control

The surfaces were then sampled using a cotton swab at the following
time-points:
![\\{0, 1, 2, 4, 7, 24\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B0%2C%201%2C%202%2C%204%2C%207%2C%2024%5C%7D "\{0, 1, 2, 4, 7, 24\}").

**NB** The control samples were diluted by $10^$3 before plating out due
to high concentrations anticipated.

![](write_up_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

------------------------------------------------------------------------

## Governing equations

We consider a number of models which potentially represent the phyisical
action of cleaning method on hands. Below we propose three such models,
with increasing complexity

-   Option 1. Re-growth on hands but no recontamination from contacts

![\\dfrac{dy}{dt}=ry-d \\exp(-gt)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdfrac%7Bdy%7D%7Bdt%7D%3Dry-d%20%5Cexp%28-gt%29 "\dfrac{dy}{dt}=ry-d \exp(-gt)")

where
![r](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r "r")
is the regrowth rate per unit time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t"),
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")
is the maximum efficacy of cleaning method and
![g](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g "g")
is the decay rate induced by the cleaning method with units
![^{-t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5E%7B-t%7D "^{-t}").

-   Option 2. Re-growth on hands and recontamination from surface
    contacts over time

![\\dfrac{dy}{dt}=l + ry-d \\exp(-gt)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdfrac%7Bdy%7D%7Bdt%7D%3Dl%20%2B%20ry-d%20%5Cexp%28-gt%29 "\dfrac{dy}{dt}=l + ry-d \exp(-gt)")
