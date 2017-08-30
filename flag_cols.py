from astropy.io import fits

def add_flag_cols(infile):
    print("reading FITS file: " , infile)
    fitsf=fits.open(infile)
    #get starting header to update and original columns
    prihdu=fitsf[0]
    toupdate_hdr=fitsf[1].header.__copy__()
    orig_cols = fitsf[1].data.columns
    #
    # Column Flag: Create a new column for a flag array to be added to the 
    # existing FITS_rec
    #
    mcols = [] 
    initarray = [0]
    mcols.append(fits.Column(name='MFLAG', format='1J', array= [0]))
    new_mcol = fits.ColDefs(mcols)
    #
    # Column Flag: Create a new column for a flag array to be added to the 
    # existing FITS_rec
    ccols = [] 
    initarray = [0]
    ccols.append(fits.Column(name='CFLAG', format='1024J', array= initarray))
    new_cols = fits.ColDefs(ccols)
    #
    #update table data with mixer flag column
    mcol_hdu = fits.BinTableHDU.from_columns(orig_cols + new_mcol)
    cols=mcol_hdu.columns
    new_hdu = fits.BinTableHDU.from_columns(cols + new_cols)
    #
    #update header with flag column information
    new_hdr=new_hdu.header.__copy__()
    #update the original header with new information of additional column
    toupdate_hdr.extend(new_hdr,strip=False,update=True)
    #add updated header to HDU
    new_hdu.header=toupdate_hdr
    #
    #make new HDU
    thdulist=fits.HDUList([fitsf[0],new_hdu])
    #
    return thdulist
    