import ROOT
import math

from ROOT import TGraphAsymmErrors

def systematics( iline, lines ):
    thislinepart = ' '.join( lines[iline].split() ).strip().split(' ')
    
    linenr = iline + 1
    nextlinepart = ' '.join( lines[linenr].split() ).strip().split(' ')
    
    uppsq = 0
    lowsq = 0
    
    while ( '+' in nextlinepart[0] ):
        uppsq += float( nextlinepart[0] )**2
        lowsq += float( nextlinepart[1] )**2
        
        linenr += 1
        nextlinepart = ' '.join( lines[linenr].split() ).strip().split(' ')
    
    return (math.sqrt( lowsq ), math.sqrt( uppsq ))


def check_number( strings ):
    for s in strings:
        try:
            float(s)
        except ValueError:
            return False
    return True

def readFile( fname , maxenergy = 4.8 ):
    
    infile  = file( fname, 'r' )
    inlines = infile.readlines()
    
    graph   = TGraphAsymmErrors()

    for iline in range(0,len(inlines)):
        line = inlines[iline]
        
        if '*' in line: continue
            
        lineparts = ' '.join( line.split() ).strip().split(' ')
            
        if not check_number( lineparts[0:5] ): continue
        
        energy  = float(lineparts[0])
        R       = float(lineparts[3])
            
        if ( energy > maxenergy ): continue
            
        statupp = abs(float(lineparts[4]))
        statlow = abs(float(lineparts[5]))
            
        (systlow, systupp) = systematics( iline, inlines )
            
            
        ipoint = graph.GetN()
                
        graph.SetPoint( ipoint, energy, R )
                
        graph.SetPointError( ipoint, 0 , 0 ,
                             math.sqrt( statlow*statlow + 1e-4*systlow*systlow*R*R ),
                             math.sqrt( statupp*statupp + 1e-4*systupp*systupp*R*R ) )

    return graph



