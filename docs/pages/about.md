# about

**_slick_** is a tool to calculate CO, [CI] and [CII] line intensities from extremely large datasets of clouds/galaxies in hydrodynamical simulations. It runs radiative transfer equations on each single particle of the simulation, and outputs:
1) A table with all clouds' and galaxies' luminosities
2) A lightcone considering a certain redshift range
3) Intensity maps for specific redshifts

It was designed to be used on [SIMBA](http://simba.roe.ac.uk/) simulation datasets, but it will be expanded and generalized in the near future.

?> We are currently implementing an option to calculate estimated luminosities even faster on extended simulation boxes using Machine Learning techniques.

### developers and collaborators

*slick* was developed by **[Karolina Garcia](https://karolinagarcia.github.io/)** in collaboration with:
- Desika Narayanan
- Gergo Popping
- Sagan Sutherland
- Thomas Greve
- Romeel Dave

And it has been in use by:
- Anirudh Ravishankar