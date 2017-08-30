#!/usr/bin/env python

import numpy as numpy
import math
import warnings

from MonteCarloError import MonteCarloError
from ConvergenceError import ConvergenceError
import Tools

__author__ = "Do Kester"
__year__ = 2016
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Do"
__status__ = "Development"

# j2py-generated source for module Fitter
#  * This file is part of Herschel Common Science System ( HCSS ).
#  * Copyright 2001-2014 Herschel Science Ground Segment Consortium
#  *
#  * HCSS is free software: you can redistribute it and/or modify
#  * it under the terms of the GNU Lesser General Public License as
#  * published by the Free Software Foundation, either version 3 of
#  * the License, or ( at your option ) any later version.
#  *
#  * HCSS is distributed in the hope that it will be useful,
#  * but WITHOUT ANY WARRANTY; without even the implied warranty of
#  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  * GNU Lesser General Public License for more details.
#  *
#  * You should have received a copy of the GNU Lesser General
#  * Public License along with HCSS.
#  * If not, see <http://www.gnu.org/licenses/>.
#  * $Id: Fitter.java,v 1.3 2014/04/09 12:17:09 do Exp $
#  *
#  *    2003 - 2014 Do Kester, SRON (JAVA source)
#  *    2016        Do Kester

class BaseFitter( object ):
    """
    Base class for all Fitters.

    The Fitter class is to be used in conjunction with *Model classes.

    The Fitter class and its descendants fit data to a model. Fitter itself
    is the variant for linear models, ie. models linear in its parameters.

    For both linear and nonlinear models it holds that once the optimal
    estimate of the parameters is found, a variety of calculations is exactly
    the same: standard deviations, noise scale, evidence and model errors.
    They all derive more or less from the inverse Hessian matrix ( aka the
    covariance matrix ). All these calculations are in this Fitter class.
    Other Fitter classes relegate their calculation in these issues to this one.

    Examples
    --------
    # It is not possible to use this class. User Fitter, CurveFitter etc. in stead

    Note Also
    ---------
    1. The calculation of the evidence is an Gaussian approximation which is
       only exact for linear models with a fixed scale.
    2. Attributes labelled as read only should not be set by a user.

    Category:    Mathematics/Fitting

    Attributes
    ----------
    model : Model
        the model to be fitted
    xdata : array_like
        independent variable(s)
    nxdata : int (read only)
        length of the xdata vector(s)
    ndim : int (read only)
        number of xdata vectors
    imastent : ImageAssistant
        to convert images to pixels indices, needed for a fit
    design : matrix (read only)
        the design matrix (partial of model to parameters)
        returns self.getDesign()
    hessian : matrix (read only)
        the hessian matrix
        returns self.getHessian()

    Attributes (available after a call to fit())
    ----------
    chisq : float (read only)
        chisquared of the fit
    stdevs, standardDeviations : array_like (read only)
        the standard deviations on the parameters
        returns self.getStandardDeviations()
    scale : float (read only)
        the noise scale
        returns self.getScale()
    covariance : matrix (read only)
        the covariance matrix
        returns self.getCovarianceMatrix()
    sumwgt : float (read only)
        sum of the weights
    logZ : float (read only)
        the e-log of the evidence
        returns self.getLogZ()
    evidence : float (read only)
        the 10log of the evidence (logZ / log(10))
        returns self.getEvidence()

    Attributes (available after a call to getLogZ() or getEvidence())
    ----------
    logOccam : float (read only)
        Occam factor
    logLikelihood : float (read only)
        log of the likelihood


    """

    #  *****CONSTRUCTORS********************************************************
    def __init__( self, xdata, model, map=False ):
        """
        Create a new Fitter, providing inputs and model.

        A Fitter class is defined by its model and the input vectors ( the
        independent variable ). When a fit to another model and/or another
        input vector is needed a new object should be created.

        Parameters
        ----------
        xdata : array_like
            independent input variable(s)
        model : Model
            the model function to be fitted
        map : bool (False)
            When true, the xdata should be interpreted as a map.
            The fitting is done on the pixel indices of the map,
            using ImageAssistant

        Raises
        ------
        ValueError when one of the following is true
            1. Dimensionality of model and input does not match.
            2. Nans in input stream.
            3. Model is not the head if the compound model chain.

        """
        if model != model._head:
            raise ValueError( "Model is not the head of a Compound model chain" )
        if  numpy.any( numpy.isnan( xdata ) ) :
            raise ValueError( "NaNs in xdata array" )

        if map :
            self.ia = ImageAssistant()
            self.xdata = self.ia.getIndices( xdata )
        else :
            self.ia = None
            self.xdata = xdata

        ninp = Tools.length( xdata )
        self.nxdata = ninp
        self.ndim = xdata.ndim
        self.model = model

        if self.ndim != model.ndim:
            raise ValueError( "Model and xdata must be of the same dimensionality." )


    #  *****FIT*****************************************************************
    def modelFit( self, ydata, weights=None ):
        """
        Return model fitted to the data.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted.
        weights : array_like
            weights to be used

        """
        self.model.parameters = self.fit( ydata, weights=weights )
        if self.model.noiseScale and not self.model.noiseScale.fixed :
            self.model.noiseScale.scale = self.scale

        yfit = self.model.result( self.xdata )
        return yfit if self.ia is None else self.ia.resizeData( yfit )

    def limitsFit( self, fitmethod, ydata, weights=None ) :
        """
        Fit the data to the model.
        When a parameter(s) transgresses the limits, it set and fixed at that limit
        and the fit is done again, excluding the parameter(s)
        When the chisq landscape is largely monomodal (no local minima) this is OK.

        Parameter
        ---------
        fitmethod : callable fitmethod( ydata, weights=weights )
            A fit method from the BaseFitter family
        ydata : array_like
            data that the model needs to be fit to
        weights : array_like
            weights partaining to the data.

        Returns
        -------
        pars : array_like
            the parameters of the fit

        Raises
        ------
        Warning when parameters have been reset at the limits.

        """
        pars = fitmethod( ydata, weights=weights )          # perform the fit
        if self.model.priors is None :                           # no priors -> no limits
            return pars

        sfix = self.model._keepFixed                        # save original fixed setting

        npchain = self.model.npchain
        params = self.model.parameters
        index = self.model.fitIndex
        if index is None :
            index = range( npchain )
            fix = []
        else :
            fix = sfix.copy()

        ool = 0
        for k in index :                                    # check limits for params
            if params[k] < self.model.getPrior( k ).lowLimit :
                fix = fix + [k]                             # add to original fixed
                params[k] = self.model.getPrior( k ).lowLimit
                ool += 1
            elif params[k] > self.model.getPrior( k ).highLimit :
                fix = fix + [k]
                params[k] = self.model.getPrior( k ).highLimit
                ool += 1

        if ool > 0 :                                        # some transgressions
            self.model.keepFixed( fix )                     # fix them
            pars = fitmethod( ydata, weights=weights )      # run fit again
            self.model.keepFixed( sfix )                    # reset original fitIndex
            if sfix is not None :                           # select parameters
                pars = [self.model.parameters[k] for k in self.model.fitIndex]
                fix = list( set( fix ) - set( sfix ) )
            warnings.warn( "Parameters ", fix, " exceeded limits." )

        return pars                                         # return parameters


    def fit( self, ydata, weights=None ) :
        """
        Return model parameters fitted to the data.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted.
        weights : array_like
            weights to be used

        Raises
        ------
        ConvergenceError. BaseFitter cannot perform fits by itself.

        """
        raise ConvergenceError( "BaseFitter is a base class, not suitable itself to perform fits." )

    def checkNan( self, ydata, weights=None ):
        if not numpy.all( numpy.isfinite( ydata ) ) :
            raise ValueError( "Fitter: NaNs or Infs in ydata" )
        if  weights is not None and not numpy.all( numpy.isfinite( weights ) ) :
            raise ValueError( "Fitter: NaNs or Infs in weights" )

    def __getattr__( self, name ) :
        """
        Return value belonging to attribute with name.

        Parameters
        ----------
        name : string
            name of the attribute
        """
        if name in ['logOccam', 'logLikelihood', 'chisq', 'sumwgt'] :
            raise AttributeError( str( self ) + ": " + name + " is not yet available." )
        elif name == 'design' :
            return self.getDesign()
        elif name == 'hessian' :
            return self.getHessian()
        elif name == 'covariance' :
            return self.getCovarianceMatrix()
        elif name == 'stdevs' or name == 'standardDeviations' :
            return self.getStandardDeviations()
        elif name == 'scale' :
            return self.getScale()
        elif name == 'evidence' :
            return self.getEvidence()
        elif name == 'logZ' :
            return self.getLogZ()
        else :
            raise AttributeError( str( self ) + ": Unknown attribute " + name )

        return None


    #  *****VECTOR**************************************************************
    def getVector( self, ydata, weights=None ):
        """
        Return the &beta;-vector.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted.
        weights : array_like
            weights to be used

        """
        dcw = ydata if weights is None else ydata * weights

        index = self.model.fitIndex
        if self.model.isNullModel() :
            return numpy.asarray( 0 )
        design = self.getDesign( index=index )
        return numpy.inner( design.transpose(), dcw )

    #  *****HESSIAN**************************************************************
    def getHessian( self, params=None, weights=None, index=None ):
        """
        Calculates the hessian matrix for a given set of model parameters.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted.
        weights : array_like
            weights to be used
        index : list of int
            index of parameters to be fixed

        """
        if params is None : params = self.model.parameters
        self.sumwgt = self.nxdata if weights is None else numpy.sum( weights )

        if self.model.isNullModel() :
            return

        design = self.getDesign( xdata=self.xdata, params=params )
        if index is not None :
            design = design[:,index]
        design = design.transpose()
        if weights is not None :
            hessian = numpy.inner( design, design * weights )
        else :
            hessian = numpy.inner( design, design )
        return hessian

#      * TBD Condition number see Wikipedia: Condition Number and Matrix Norm

    #  *************************************************************************
    def getInverseHessian( self, params=None, weights=None, index=None ):
        """
        Return the inverse of the Hessian Matrix, H.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted.
        weights : array_like
            weights to be used
        index : list of int
            index of parameters to be fixed

        """
        return numpy.linalg.inv( self.getHessian( params, weights, index ) )

    def getCovarianceMatrix( self ):
        """
        Returns the inverse hessian matrix over the fitted parameters,
                multiplied by the variance.

        Stdevs are found from this as np.sqrt( np.diag( covarianceMatrix ) )

        """
        return self.getInverseHessian( index=self.model.fitIndex ) * self.makeVariance()

    #  *****DESIGN**************************************************************
    def getDesign( self, params=None, xdata=None, index=None ):
        """
        Return the design matrix, D.
        The design matrix is also known as the Jacobian Matrix.

        Parameters
        ----------
        xdata : array_like
            the independent input data
        params : array_like
            parameters of the model
        index : list of int
            index of parameters to be fixed

        """
        if params is None : params = self.model.parameters
        if xdata is None :  xdata = self.xdata

        design = self.model.partial( xdata, params )
        if index is not None :
            design = design[:,index]
        return design

    #  *****CHI-SQUARED*********************************************************
    def chiSquared( self, ydata, params=None, weights=None ):
        """
        Calculates Chi-Squared for data and weights.

        Parameters
        ----------
        ydata : array_like
            the data vector to be fitted.
        params : array_like
            parameters for the model
        weights : array_like
            weights to be used

        """

        res2 = numpy.square( ydata - self.model.result( self.xdata, params ) )
        if weights is not None:
            res2 *= weights
            self.sumwgt = numpy.sum( weights )
        else:
            self.sumwgt = self.nxdata
        self.chisq = numpy.sum( res2 )
        if self.chisq <= 0 :
            raise ValueError( str( self ) + ": chisq <= 0" )
        return self.chisq

    #  *****STANDARD DEVIATION**************************************************
    def getStandardDeviations( self ):
        """
        Calculates of standard deviations pertaining to the parameters.

        .. math::
            \sigma_i =  * \sqrt( C_{i,i} )

        where C is the Covariance matrix, the inverse of the Hessian Matrix and
        s is the noiseScale.

        Standard deviation are calculated for the fitted parameters only.

        Note that the stdev will decrease with sqrt( N ) of the number of
        datapoints while the noise scale, s, does not.

        """
        nfit = self.model.getNumberOfFittedParameters( )
        matrix = self.getInverseHessian( index=self.model.fitIndex )
        self.model.stdevs = numpy.sqrt( matrix.diagonal() * self.chisq / ( self.nxdata - nfit ) )

        if self.model.noiseScale and not self.model.noiseScale.fixed :
            self.model.noiseScale.stdevScale = self.scale / math.sqrt( 2 * self.nxdata )

        return self.model.stdevs

    #  *****MONTE CARLO ERROR***************************************************
    def monteCarloError( self, xdata=None, monteCarlo=None):
        """
        Calculates :math:\sigma:math:-confidence regions on the model given some inputs.

        From the full covariance matrix (inverse of the Hessian) random
        samples are drawn, which are added to the parameters. With this new
        set of parameters the model is calculated. This procedure is done
        by default, 25 times.
        The standard deviation of the models is returned as the error bar.

        The calculation of the confidence region is delegated to the class
        MonteCarloError. For tweaking of that class can be done outside BaseFitter.

        Parameters
        ----------
        xdata : array_like
            input data over which to calculate the error bars.
        monteCarlo : MonteCarloError
            a ready-made MonteCarloError class.

        """
        if xdata is None : xdata = self.xdata
        if monteCarlo is None :
            monteCarlo = MonteCarloError( xdata, self.model, self.chisq,
                                          covariance=self.getCovarianceMatrix( ) )
        return monteCarlo.getError( xdata )

    #  *************************************************************************
    def makeVariance( self ):
        dof = self.nxdata - self.model.getNumberOfFittedParameters()
        scale = 1.0 if self.model.noiseScale is None else self.model.noiseScale.scale
        chisq = self.chisq
        if self.minimumScale() :
            chisq += self.sumwgt * scale * scale
        var = chisq / dof if self.autoScale() else scale * scale
        return var

    def autoScale( self ) :
        return self.model.autoScale()

    def minimumScale( self ) :
        return self.model.minimumScale()


    #  *****EVIDENCE************************************************************
    def getEvidence( self, limits=None, noiseLimits=None ):
        """
        Calculation of the evidence, log10( Z ), for the model given the data.

        ../math::
            E = \log10( P( Model | data ) )

        The calculation of the evidence uses a Gaussion approximation of the Posterior
        probability.
        It needs to know the limits of the parameters (and the noise scale if applicable),
        either from the priors in the model or from keywords "limits/noiseLimits".


        Parameters
        ----------
        limits : list of 2 floats/array_likes
            possible range of the parameters. ( [low,high] )
        noiseLimits : list of 2 floats
            possible range on noise scale ( [low,high] )

        Raises
        ------
        ValueError when no Prior is available

        """
        return self.getLogZ( limits, noiseLimits ) / math.log( 10.0 )

    def getLogZ( self, limits=None, noiseLimits=None ):
        """
        Calculation of the evidence, log( Z ), for the model given the data.

        .. math::
            logZ = \log( P( Model | data ) )

        The calculation of the evidence uses a Gaussion approximation of the Posterior
        probability.
        It needs to know the limits of the parameters (and the noise scale if applicable),
        either from the priors in the model or from keywords "limits/noiseLimits".


        Parameters
        ----------
        limits : list of 2 floats/array_likes
            possible range of the parameters. ( [low,high] )
        noiseLimits : list of 2 floats
            possible range on noise scale ( [low,high] )

        Raises
        ------
        ValueError when no Prior is available

        """
        hasNoPrior = 0
        priors = self.model.priors
        if limits is not None :
            prirange = limits[1] - limits[0]
        elif self.model.hasLimits() :
            prirange = [p.getIntegral() for p in priors]        # convert to list
        else :
            hasNoPrior += 1

        if noiseLimits is not None :
            scalerange = math.log( noiseLimits[1] ) - math.log( noiseLimits[0] )
        elif self.model.hasNoiseLimits() :
            scalerange = self.model.noiseScale.prior.getIntegral()
        elif self.autoScale() :
            hasNoPrior += 2

        if hasNoPrior > 0 :
            if hasNoPrior == 1 : mess = "parameters."
            if hasNoPrior == 2 : mess = "noise scale."
            if hasNoPrior == 3 : mess = "parameters and noise scale."
            raise ValueError( "No limits provided on " + mess + " Cannot calculate evidence." )

        npar = self.model.getNumberOfParameters( )
        nfit = self.model.getNumberOfFittedParameters( )

        priorlength = Tools.length( prirange )                         # maybe less than npar
        prirange = numpy.log( numpy.asarray( prirange ) )

        if self.model.fitIndex is None :
            spr = prirange.sum()
        else :
            fi = self.model.fitIndex
            q = fi.where( fi < priorlength )
            fq = fi[q]
            priorlength = len( fq )
            spr = numpy.sum( prirange[fq] )

        if priorlength < npar :                             # add enough times the last one
            spr += ( npar - priorlength ) * ( prirange if priorlength == 1 else prirange[-1] )

        # add integral over the scale prior, if scale is to be fitted.
        if not ( self.model.noiseScale and self.model.noiseScale.fixed ) :
            spr +=  math.log( scalerange )

        self.logOccam = 0.0
        s2 = self.makeVariance( )
        if not self.model.isNullModel( ):
            hes = self.getHessian( index=self.model.fitIndex )
            logdet = math.log( numpy.linalg.det( hes ) )
            self.logOccam = -spr + 0.5 * ( nfit * math.log( 2 * math.pi * s2 ) - logdet )

        self.logLikelihood = -0.5 * ( self.sumwgt * math.log( 2 * math.pi * s2 ) + self.chisq / s2 )

        return self.logLikelihood + self.logOccam

    #  *************************************************************************
    def getScale( self ):
        """
        Return the noise scale: sqrt( chisq / DOF ).

        """
        dof = self.nxdata - self.model.getNumberOfFittedParameters()
        return math.sqrt( self.chisq / dof )

    def __str__( self ):
        """ Return name of the fitter.  """
        return "BaseFitter"

