import ROOT
import math

from math import pow

from lhcbStyle import lhcbStyle
lhcbStyle()

from ROOT import TMath

_alpha_mZ       = 1./128.
_sin_sq_theta_W = 0.23
_gz_ge          = 1./( math.sqrt( _sin_sq_theta_W*( 1. - _sin_sq_theta_W ) ) )
_mZ             = 91.2
_gZ             = 2.5


'''
    Definitions: 
    
    mW^2 = mZ^2 cos^2 theta_W
    
    GF = pi*alpha/[ sqrt(2) . mW^2 . sin^2 theta_W ] 
       = pi*alpha/[ sqrt(2) . mZ^2 cos^2 theta_W . sin^2 theta_W ]
       
    GF . mZ^2 / [ 8 . pi . sqrt(2)] = alpha/[ 16 . cos^2 theta_W . sin^2 theta_W ]
'''



class fermion:
    '''
        Simple class holding fermion properties 
    '''
    def __init__( self,
                  fermion_name,
                  fermion_charge,
                  fermion_weak_isospin,
                  fermion_mass ):
        self.name = fermion_name
        self.q    = fermion_charge
        self.m    = fermion_mass
        self.t3   = fermion_weak_isospin

    def gv( self ):
        global _sin_sq_theta_W
        return self.t3 - 2.*self.q*_sin_sq_theta_W
    
    def ga( self ):
        return self.t3
    
    def active( self, s ):
        return ( math.sqrt(s) > 2.*self.m )

    def print_properties( self ):
        print self.name
        print  ( ' q = ' + str(self.q) ,
                ' t3 = ' + str(self.t3) ,
                ' gA = ' + str(self.ga()) ,
                ' gV = ' + str(self.gv()) )



def alpha( s, alpha_mZ, mZ ):
    '''
        Compute running of alpha
    '''
    return alpha_mZ/(1. - ( alpha_mZ/(3.*TMath.Pi()) )*math.log( s/ mZ*mZ ) )


def chi_real( s, mZ, gZ, sin_sq_theta_W ):
    '''
        Evaluate Re(chi)
    '''
    result = ( ( s*s - s*mZ*mZ )/( math.pow( s - mZ*mZ , 2 ) + mZ*mZ*gZ*gZ ) )
    return result/( 16.*sin_sq_theta_W*( 1. - sin_sq_theta_W ) )

def chi_square( s, mZ, gZ, sin_sq_theta_W ):
    '''
        Evaluate chi^2
    '''
    result = ( s*s )/( math.pow( s - mZ*mZ , 2 ) + mZ*mZ*gZ*gZ )
    return result/pow( 16.*sin_sq_theta_W*( 1. - sin_sq_theta_W ) , 2 )


def differential_decay_rate( s, i, o, NC = 1. ):
    '''
        Differential decay rate as function of s with gamma-Z
    '''

    global _mZ
    global _gZ
    global _alpha_mZ
    global _sin_sq_theta_W
    
    result = 0
    if not i.active( s ):
        return 0

    for p in o:
        if not p.active( s ): continue
        
        result += ( pow( i.q*p.q, 2 ) )
        
        result += ( 8.*i.q*p.q*i.gv()*p.gv()*chi_real(s, _mZ, _gZ, _sin_sq_theta_W ) )
        
        result += ( 16.*( pow(i.ga(),2) + pow(i.gv(),2) )*( pow(p.ga(),2) + pow(p.gv(),2) )*
                    chi_square(s, _mZ, _gZ, _sin_sq_theta_W) )
    
    return (8./3.)*NC*TMath.Pi()*result*pow( alpha(s,_alpha_mZ,_mZ), 2 )/(2.*s)


def differential_decay_rate_qed( s, NC = 1. ):
    '''
        Basic QED differential decay rate
    '''
    return 4.*NC*TMath.Pi()*pow( alpha(s,_alpha_mZ,_mZ), 2 )/(3.*s)

def gamma_Z_interference( s, i, o, NC = 1 ):
    '''
        Interference term between g and Z0
    '''
    global _mZ
    global _gZ
    global _alpha_mZ
    global _sin_sq_theta_W
    
    result = 0
    if not i.active( s ):
        return 0

    for p in o:
        if not p.active( s ): continue
        
        result += ( 8.*i.q*p.q*i.gv()*p.gv()*chi_real(s, _mZ, _gZ, _sin_sq_theta_W ) )

    return (8./3.)*NC*TMath.Pi()*result*pow( alpha(s,_alpha_mZ,_mZ), 2 )/(2.*s)



def forward_backward_asymmetry( s, i, o ):
    '''
        Forward-backward asymmetry
    '''
    global _mZ
    global _gZ
    global _alpha_mZ
    global _sin_sq_theta_W
    
    if i.active( s ) and o.active( s ):
        result = 0
        result += ( 16.*i.q*o.q*i.ga()*o.ga()*chi_real( s, _mZ, _gZ, _sin_sq_theta_W ) )
        result += ( 16.*8.*i.gv()*o.gv()*i.ga()*o.ga()*chi_square( s, _mZ, _gZ, _sin_sq_theta_W ) )
    
        result = result*pow( alpha( s, _alpha_mZ, _mZ ), 2 )*TMath.Pi()/(2.*s)
    
        return result/differential_decay_rate( s, i, [o] )
    return 0


def afb_with_sqrts( i, o ):
    '''
    Create graph showing the AFB
    '''
    sqrts = 0.2
    point = 0

    result = TGraph()

    
    while ( sqrts < 160. ):
        afb = forward_backward_asymmetry( sqrts**2, i, o )
        result.SetPoint( point, sqrts, afb )

        
        sqrts += 0.1
        point += 1

    
    return result


def interference_with_sqrts( i, o, NC = 1 ):
    '''
        Create graph showing interference
    '''
    sqrts = 0.2
    point = 0

    
    result = TGraph()

    
    while ( sqrts < 160. ):
        ratio = gamma_Z_interference( sqrts**2, i, o, NC )
        result.SetPoint( point, sqrts, ratio/differential_decay_rate_qed( sqrts**2 ) )

        
        sqrts += 0.1
        point += 1

    
    return result




def ratio_with_sqrts( i, o, NC = 1. ):
    '''
        Create graph showing the R ratio
    '''
    sqrts = 0.2
    point = 0

    result = TGraph()

    while ( sqrts < 160. ):
        ratio= differential_decay_rate( sqrts**2, i, o, NC )/differential_decay_rate_qed( sqrts**2 )
        
        result.SetPoint( point, sqrts, ratio )
    
        sqrts += 0.1
        point += 1
    
    return result





electron = fermion( 'electron' , -1, -1./2., 5.11e-4 )
muon     = fermion( 'muon'     , -1, -1./2., 0.106 )

quarks   = [
    fermion( 'up'     ,  2./3. ,  1./2., 2.3e-3 ),
    fermion( 'down'   , -1./3. , -1./2., 4.8e-3 ),
    fermion( 'charm'  ,  2./3. ,  1./2., 1.3 ),
    fermion( 'strange', -1./3. , -1./2., 9.5e-2 ),
    fermion( 'top'    ,  2./3. ,  1./2., 173. ),
    fermion( 'bottom' , -1./3. , -1./2., 4.7 )
    ]

from ROOT import TGraph

graph_m = TGraph()
graph_m_qed = TGraph()
graph_i = TGraph()
graph_h = TGraph()
graph_r = TGraph()


graph_afb = TGraph()

from pyDatasetForRratio import readFile

graph_data = readFile( fname = 'rpp2012-hadronicrpp_page1001.dat.txt', maxenergy = 160.0 )

sqrts = 0.2
point = 0

while ( sqrts < 160.0 ):
    
    rate_m_full = differential_decay_rate( sqrts**2 , electron , [muon] )
    rate_m_qed  = differential_decay_rate_qed( sqrts**2 )
    
    rate_h      = differential_decay_rate( sqrts**2 , electron , quarks, NC = 3. )

    rate_i      = gamma_Z_interference( sqrts**2, electron, quarks, NC = 3. )
    
    afb         = forward_backward_asymmetry( sqrts**2, electron, muon )
    
    graph_afb.SetPoint( graph_afb.GetN(), sqrts, afb )
    
    if ( rate_h > 0 and rate_m_qed > 0 ):
        
        graph_m.SetPoint( point, sqrts, rate_m_full )
        graph_m_qed.SetPoint( point, sqrts, rate_m_qed )
        graph_h.SetPoint( point, sqrts, rate_h )
        graph_r.SetPoint( point, sqrts, rate_h/rate_m_qed )
        graph_i.SetPoint( point, sqrts, rate_i/rate_m_qed )
    
    sqrts += 0.1
    point += 1


from ROOT import gPad, TCanvas

can_r = TCanvas('can_r')

graph_r.Draw('AL')
graph_r.SetLineColor( ROOT.kRed )
graph_r.GetXaxis().SetTitle( '#sqrt{s} [GeV]' )
graph_r.GetYaxis().SetTitle( '#it{R}' )

gPad.SetLogy()

graph_data.Draw('P+')

can_afb = TCanvas('can_afb')
graph_afb.Draw('AL')

graph_afb.GetXaxis().SetTitle( '#sqrt{s} [GeV]' )
graph_afb.GetYaxis().SetTitle( '#it{A}_{FB}' )



'''
    Change value of sin^2 theta_W
'''

'''
v23 = ratio_with_sqrts( electron, quarks, NC = 3. )

_sin_sq_theta_W = 0.25
v25 = ratio_with_sqrts( electron, quarks, NC = 3. )

from ROOT import TFile
tfile = TFile( 'R.root' , 'RECREATE' )
tfile.cd() 

v23.Write('0.23')
v25.Write('0.25')
tfile.Close()
'''

