from astropy.io import fits

def addRowFlag(infile,flagcolname,initvalue):
    #
    #This task will add a row column to an existing FITS table (infile) and call that column "flaghdrname"
    #
    print("reading FITS file: " , infile)
    fitsf=fits.open(infile)
    #get starting header to update and original columns
    prihdu=fitsf[0]
    #
    # Get a copy of the FITS header
    #
    toupdate_hdr=fitsf[1].header.__copy__()
    #
    # Get the original list of columns
    #
    orig_cols = fitsf[1].data.columns
    #
    # Row Flag: Create a new column for an integration/row flag to be added to the 
    # existing FITS_rec.
    #
    # The name of new column is "flagcolname"
    #
    # create an empty row flag 
    mcols = [] 
    #
    #instantiate to 0
    #
    initarray=[0]
    initarray[0] = initvalue 
    mcols.append(fits.Column(name=flagcolname, format='1J', array= initarray))
    new_mcol = fits.ColDefs(mcols)
    #
    #update table data with mixer flag column
    mcol_hdu = fits.BinTableHDU.from_columns(orig_cols + new_mcol)
    cols=mcol_hdu.columns
    new_hdu = fits.BinTableHDU.from_columns(cols + new_cols)
    #
    #make new HDU
    thdulist=fits.HDUList([fitsf[0],new_hdu])
    #
    return thdulist
    
