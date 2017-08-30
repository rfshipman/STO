#this script repojects an RA/DEC image onto galactic coordinates.  
spfile='spire500.fits'
hdu1=fits.open(spfile)[1]
hdu2=fits.open("skv9919108574256.fits")[0]


from astropy.wcs import WCS
import matplotlib.pyplot as plt

ax1 = plt.subplot(1,2,1, projection=WCS(hdu1.header))
ax1.imshow(hdu1.data, origin='lower', vmin=0, vmax=10.)
ax1.coords.grid(color='white')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')
ax1.set_title('SPIRE 500 micron')

ax2 = plt.subplot(1,2,2, projection=WCS(hdu2.header))
ax2.imshow(hdu2.data, origin='lower', vmin=200, vmax=2000)
ax2.coords.grid(color='white')
ax2.coords['glon'].set_axislabel('Galactic Longitude')
ax2.coords['glat'].set_axislabel('Galactic Latitude')
ax2.coords['glat'].set_axislabel_position('r')
ax2.coords['glat'].set_ticklabel_position('r')
ax2.set_title('IRAS 100')


from reproject import reproject_interp
array, footprint = reproject_interp(hdu1, hdu2.header)

#from astropy.wcs import WCS
#import matplotlib.pyplot as plt

ax1 = plt.subplot(1,2,1, projection=WCS(hdu2.header))
ax1.imshow(array, origin='lower', vmin=0.0, vmax=15)
ax1.coords.grid(color='white')
ax1.coords['glon'].set_axislabel('Galactic lon')
ax1.coords['glat'].set_axislabel('Galactic lat')
ax1.set_title('Reprojected SPIRE image')

ax2 = plt.subplot(1,2,2, projection=WCS(hdu2.header))
ax2.imshow(footprint, origin='lower', vmin=-10, vmax=1.5)
ax2.coords.grid(color='white')
ax1.coords['glon'].set_axislabel('Galactic lon')
ax1.coords['glat'].set_axislabel('Galactic lat')
ax2.coords['glat'].set_axislabel_position('r')
ax2.coords['glat'].set_ticklabel_position('r')
ax2.set_title('SPIRE footprint')

fits.writeto('spire_galctic_header.fits', array, hdu2.header, clobber=True)