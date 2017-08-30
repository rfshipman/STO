import Fitter
import SineModel
import PolynomialModel 
import SineAmpModel
import numpy as np
from pylab import *
from scipy import fftpack
import unittest
import logging



def find_fringes(x, y, w):
    logger = logging.getLogger('Fitter.Fitter')
    #x = x data
    #y = y dependent variable
    #w = weights (0-1) 1 = use, 0 = don't use
    #first time through only assume a constant polynomial model
    model=PolynomialModel.PolynomialModel(0)
    testmodel=PolynomialModel.PolynomialModel(0)
    #use the maximum of the data as limit
    ymax=max(y)

    model.setLimits(lowLimits=[-ymax], highLimits=[ymax])
    model.setNoiseLimits(lowLimit=0.00001, highLimit=ymax)
    testmodel.setLimits(lowLimits=[-ymax], highLimits=[ymax])
    testmodel.setNoiseLimits(lowLimit=0.00001, highLimit=ymax)
    # fit the intital offset , just to get a baseline evidence measure
    fitter=Fitter.Fitter(x,model)
    testfitter=Fitter.Fitter(x,testmodel)
    params=testfitter.fit(y,w)
    lastevidence=testfitter.evidence
    testevidence=lastevidence
    #
    # create and keep an original data vector
    yori=y.copy()
    #
    print(testevidence)
    #remove polynomial model from the data
    y = yori - model(x)
    # create some variables to keep track of fringes
    stw_d=[]
    #
    try:
        while testevidence >= lastevidence:
            p=find_powpeak(x,y)
            sinetestmodel=SineAmpModel.SineAmpModel(p)
            sinetestmodel.setLimits(lowLimits=[-ymax,-ymax], highLimits=[ymax,ymax])
            sinetestmodel.setNoiseLimits(lowLimit=0.00001, highLimit=ymax)
            testmodel.addModel(sinetestmodel) 
            testfitter=Fitter.Fitter(x,testmodel)
            testparams=testfitter.fit(yori,w)
            testevidence=testfitter.evidence
            print("test:  ", testevidence,"previous:  ", lastevidence)
            print(testevidence > lastevidence)
            if testevidence > lastevidence:
                sinemodel=SineAmpModel.SineAmpModel(p)
                sinemodel.setLimits(lowLimits=[-ymax,-ymax], highLimits=[ymax,ymax])
                sinemodel.setNoiseLimits(lowLimit=0.00001, highLimit=ymax)
                #print("test model:  ",testmodel)
                #
                model.addModel(sinemodel)
                #print("Sine alone:  ",sinemodel)
                #print("OK model:  ",model)
                fitter=Fitter.Fitter(x,model)
                params=fitter.fit(yori,w)
                #Update the evidence with a good model
                lastevidence=fitter.evidence
                #subtract the good model from the data
                y=yori-model(x)
                #p1=plt.plot(x,yori,color="Blue")
                #p2=plt.plot(x,model(x),color="Red")
                #plt.show()
                distance=300.0/2.0/(1000.0/p)
                stw_d.append(distance)
                #testmodel=model
    except Exception as e:
        logger.error("there was an exception %s", str(e))
    print("Done")
    return model,stw_d

    

    
def find_powpeak(x, y):
    #print(y.shape[-1])
    dx=x[1]-x[0]
    #
    #Find an initial guess at the number of cycles using a power spectrum
    ft=fftpack.rfft(y)
    kfft=fftpack.rfftfreq(x.shape[-1],dx)
    #don't use the zero frequency
    k=kfft[1:]
    aft=abs(ft*ft.conjugate())
    qp=np.where(aft[1:] >= max(aft[1:]))
    kpeak=float(k[qp][0].__abs__())
    #print("k peak:  ",kpeak)
    p0init=float(kpeak)
    if p0init != 0.0:
        p0range=np.arange(0.9*p0init,1.1*p0init,0.03)
        chis=[]
        #test each sine model over a range of p0 for lowest Chi2
        #print(p0range)
        for p0 in p0range:
            testmod = SineAmpModel.SineAmpModel(p0)
            linefit = Fitter.Fitter(x,testmod)
            testparam=linefit.fit(y)
            chis.append(linefit.chisq)
        #lowest Chisquared
        p0=p0range[chis.index(min(chis))]
    else:
        p0=0.0
    return p0
