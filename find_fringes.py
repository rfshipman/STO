import Fitter
import SineModel
import PolynomialModel 
import SineAmpModel
import numpy as np
from pylab import *
from scipy import signal

#import /define find_peakpow can only do this with an actual module.
#%run Find_stw.py

def find_fringes(x,y,w):
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
                p1=plt.plot(x,yori,color="Blue")
                p2=plt.plot(x,model(x),color="Red")
                plt.show()
                distance=300.0/2.0/(1000.0/p)
                stw_d.append(distance)
                #testmodel=model
    print("Done")
    return model,stw_d
