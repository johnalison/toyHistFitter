import ROOT
ROOT.gROOT.SetBatch(True)
import ROOTHelp.FancyROOTStyle
import sys, os
import numpy as np

from ROOTHelp.Utils import drawStackCompRatio

import optparse
parser = optparse.OptionParser()
parser.add_option('-o', '--outputDir',           help="OutputDir")
o, a = parser.parse_args()

if not os.path.isdir(o.outputDir):
    os.mkdir(o.outputDir)

path = o.outputDir+"/"

#
# Pointers to files
#

signalFileName = "ZH4b_nEvts_1932400_hists_nonMixed.root"
bkgFileName = "data18_wSvB_2018_hists_nonMixed.root"
bkgFileMixedName = "4b_Mixed4b_wSvB_2018_hists.root"

signalInjectedInfo = [
    {"trueFrac" :0.423614150981,
     "bkgModel":"data18_4bEmulatedwMCBranches_ZH4b_nEvts_1932400_histsMixed_3bAMix4b.root",
     "psData"  :"data18_4bEmulatedwMCBranches_ZH4b_nEvts_1932400_hists_nonMixed.root"
     },
    {"trueFrac" :0.646274452341,
     "bkgModel":"data18_4bEmulatedwMCBranches_ZH4b_nEvts_4831000_histsMixed_3bAMix4b.root",
     "psData"  :"data18_4bEmulatedwMCBranches_ZH4b_nEvts_4831000_hists_nonMixed.root"
     },
    ]





rand = ROOT.TRandom3()
ROOT.gROOT.SetBatch(True)

def getHist(inFile,histName="passXWt/threeTag/mainView/ZHSR/v4j/m_l"):
    h = inFile.Get(histName).Clone()
    return h


def signalSFToFrac(sf,hSig,hBkg):
    nSig = hSig.Integral()
    nBkg = hBkg.Integral()
    return sf*nSig / (sf*nSig+nBkg)

def makeFakeData(bkg,sig,SF):


    hFakeData = bkg.Clone()

    hFakeData.Add(sig,SF)
    hFakeData = fluctuate(hFakeData)
    hFakeData.SetName("FakeData")
    hFakeData.SetTitle("FakeData")
    return hFakeData




def fluctuate(hist):
    hist2 = hist.Clone()
    for i in range(hist.GetNbinsX()):
        mean = hist.GetBinContent(i)
        err = np.sqrt(mean)
        val = rand.Gaus(np.double(mean),np.double(err))
        while val < 0:
            val = rand.Gaus(np.double(mean),np.double(err))
        hist2.SetBinContent(i,val)
    return hist2

#print(tf1_sPlusB_func.Eval(5))



fBkg = ROOT.TFile(bkgFileName,"READ")
fSig = ROOT.TFile(signalFileName,"READ")
fBkgMixed = ROOT.TFile(bkgFileMixedName,"READ")


hBkg      = getHist(fBkg,"passXWt/threeTag/mainView/ZHSR/v4j/m_l")
hBkgMixed = getHist(fBkgMixed,"passXWt/fourTag/mainView/ZHSR/v4j/m_l")
hSig      = getHist(fSig,"passXWt/fourTag/mainView/ZHSR/v4j/m_l")


S = hSig.Integral()
B = hBkg.Integral()


can = ROOT.TCanvas("canSB","canSB")
can.cd()
hBkg.Draw()
hSig.Draw("same")
can.SaveAs(path+"RawSandB.pdf")


can.cd()
hSig.SetMarkerStyle(ROOT.kDot)
hSig.SetMarkerColor(ROOT.kBlue)
hSig.SetTitle("Signal")
hSig.SetStats(ROOT.kFALSE)
hSig.Draw()
can.SaveAs(path+"Signal.pdf")

can.cd()
hBkg.SetMarkerStyle(ROOT.kDot)
hBkg.SetMarkerColor(ROOT.kBlue)
hBkg.SetTitle("Background")
hBkg.SetStats(ROOT.kFALSE)
hBkg.Draw()
can.SaveAs(path+"Background.pdf")





#fracFit.SetCanExtend(ROOT.TH1.kAllAxes)


def doInjectionFits(name, hDataInput, hSigModel, hBkgModel, nToys):

    def sPlusB_func(x,par):    
        nBkg = hBkgModel.GetBinContent(hBkgModel.FindBin(x[0]))
        nSig = hSigModel.GetBinContent(hSigModel.FindBin(x[0]))
        return nBkg+par[0]*nSig


    tf1_sPlusB_func = ROOT.TF1("tf1_sPlusB",sPlusB_func,-10,10,1)
    tf1_sPlusB_func.SetParameter(0,1)


    sfFitted = ROOT.TH1F("sffit"+name,"Fitted vs",100,-1,2)

    #
    # Toys
    #
    for j in range(100): 
        hData = fluctuate(hDataInput)

        tf1_sPlusB_func.SetParameter(0,0)

        if j == 0:
            can.cd()
            hData.Draw()
            can.SaveAs(path+"FakeData"+str(name)+".pdf")


        fakefit = hData.Fit(tf1_sPlusB_func,"Q")
        if j == 0:
            can.SaveAs(path+"DataWithFit"+str(name)+".pdf") 
            
        sfFit    = tf1_sPlusB_func.GetParameter(0)
        sfFitErr = tf1_sPlusB_func.GetParError(0)
        #print signalSF,": ", sfFit,"+/-",sfFitErr

        sfFitted.Fill(sfFit)

    return sfFitted



def doInjectionStudy(name, hTrueBkg, hSig, hBkgModel, SFs = [0,0.02,0.03,0.1,0.2,0.4,0.6,1]):

    gr = ROOT.TGraphErrors(len(SFs))
    grBias = ROOT.TGraphErrors(len(SFs))

    fracsFit = []
    

    #
    #
    #
    for iSF, signalSF in enumerate(SFs):

        sfName = name+"_"+str(signalSF)

        #fracsFit.append(ROOT.TH1F("fracfit"+sfName,"Fitted vs",50,-1,1))
    
        hDataRaw = hTrueBkg.Clone()
        hDataRaw.Add(hSig,signalSF)
    
        #nSig = hSig.Integral()


        sFrac = signalSFToFrac(signalSF, hSig, hTrueBkg)

        #
        #  Closure test
        #
        nData = hDataRaw.Integral()
        nTrueBkg = hTrueBkg.Integral()
        sFracCalc = (nData - nTrueBkg) / nData
        print "true sFrac = ",sFrac, "(",sFracCalc,")"

        fracsFit.append(doInjectionFits(name=sfName,hDataInput=hDataRaw, hSigModel=hSig, hBkgModel=hBkgModel, nToys=100))
    
        can.cd()
        fracsFit[-1].Draw()
        can.SaveAs(path+"FracFit"+sfName+".pdf")
        
        #
        # Calculate mean fitted fraction
        #
        meanSFFit = fracsFit[-1].GetMean()
        sFracFit = signalSFToFrac(meanSFFit, hSig, hBkgModel)

        gr.SetPoint(iSF, sFracCalc, sFracFit)
        gr.SetPointError(iSF, 0, 0.1*sFracFit)

        if sFracCalc:
            sFracBias = (sFracFit-sFracCalc)/sFracCalc
        else:
            sFracBias = 0
        grBias.SetPoint(iSF, sFracCalc, sFracBias)
        grBias.SetPointError(iSF, 0, 0.1*sFracBias)

        #
        #  Draw the fit
        #
        can.cd()
        hSigFit = hSig.Clone()
        hSigFit.Scale(meanSFFit)

        drawStackCompRatio("Fit"+sfName,(hDataRaw,"PS Data"),
                           [(hBkgModel,"Background",ROOT.kAzure-9),
                            (hSigFit,"Signal Fit",ROOT.kYellow)]
                           ,yTitle="Entries",xTitle=hDataRaw.GetXaxis().GetTitle(),rTitle="Data/Fit",setLogy=0,outDir=path,cmsText="Work in Progress",lumiText="")


        print "\t fitted sFrac = ",sFracFit, "( True sFrac =",sFracCalc,")"
    return gr, grBias

#
# First do closure Test
#
grNom, grNomBias = doInjectionStudy("", hTrueBkg = hBkg, hSig = hSig, hBkgModel = hBkg) 

#
# Now closure Test with mixed bkg as model
#
grMixed, grMixedBias = doInjectionStudy("_mixed", hTrueBkg = hBkg, hSig = hSig, hBkgModel = hBkgMixed) 


grSigInj = ROOT.TGraphErrors(len(signalInjectedInfo))
grSigInjBias = ROOT.TGraphErrors(len(signalInjectedInfo))
for iter_sInfo, sInfo in enumerate(signalInjectedInfo):

    trueFracName = str(round(sInfo["trueFrac"],2))
    
    bkgModel = ROOT.TFile(sInfo["bkgModel"],"READ")
    psData   = ROOT.TFile(sInfo["psData"],"READ")

    hBkgModel = getHist(bkgModel,"passXWt/fourTag/mainView/ZHSR/v4j/m_l")
    hPSData   = getHist(psData,"passXWt/fourTag/mainView/ZHSR/v4j/m_l")

    NPSData = hPSData.Integral()
    NSig = NPSData - B 
    fracTrueCal = NSig/NPSData

    print "Truth",fracTrueCal, sInfo["trueFrac"]

    sfFitted = doInjectionFits(name="mixed_sf"+trueFracName,hDataInput=hPSData, hSigModel=hSig, hBkgModel=hBkgModel, nToys=100)
    sfFittedMean = sfFitted.GetMean()
    
    fracFitted = signalSFToFrac(sfFittedMean, hSig, hBkgModel)
    print sInfo["trueFrac"], "vs", fracFitted

    grSigInj.SetPoint(iter_sInfo, fracTrueCal, fracFitted)
    grSigInj.SetPointError(iter_sInfo, 0, 0.1*fracFitted)

    if fracTrueCal:
        sFracBias = (fracFitted-fracTrueCal)/fracTrueCal
    else:
        sFracBias = 0

    grSigInjBias.SetPoint(iter_sInfo, fracTrueCal, sFracBias)
    grSigInjBias.SetPointError(iter_sInfo, 0, 0.1*sFracBias)


    can.cd()
    sfFitted.Draw()
    can.SaveAs(path+"FracFittedmixed_sf"+trueFracName+".pdf")


    #
    #  Draw the fit
    #
    can.cd()
    hSigFit = hSig.Clone()
    hSigFit.Scale(sfFittedMean)

    drawStackCompRatio("Fit_sf"+trueFracName,(hPSData,"PS Data"),
                       [(hBkgModel,"Background",ROOT.kAzure-9),
                        (hSigFit,"Signal Fit",ROOT.kYellow)]
                       ,yTitle="Entries",xTitle=hPSData.GetXaxis().GetTitle(),rTitle="Data/Fit",setLogy=0,outDir=path,cmsText="Work in Progress",lumiText="")
    
    



hAxis = ROOT.TH1F("axis","axis;Injected Signal Fraction;Fitted Signal Fraction",10,0,0.75)
hAxis.GetYaxis().SetRangeUser(0,0.75)

lineAtOne = ROOT.TF1("tf1_sPlusB","x",-10,10,1)
lineAtOne.SetLineStyle(ROOT.kDashed)

can.cd()
hAxis.Draw()
lineAtOne.Draw("same")
grNom.Draw("P")
grMixed.SetLineColor(ROOT.kBlue)
grMixed.SetMarkerColor(ROOT.kBlue)
grMixed.Draw("P")
grSigInj.SetMarkerColor(ROOT.kRed)
grSigInj.Draw("P")

can.SaveAs(path+"fitVsInjected.pdf")



hAxisBias = ROOT.TH1F("axisBais","axisBais;Injected Signal Fraction;Signal Bias (%)",10,0,0.75)
hAxisBias.GetYaxis().SetRangeUser(-1,1)

lineAtZero = ROOT.TF1("tf1_sPlusB","0",-10,10,1)
lineAtZero.SetLineStyle(ROOT.kDashed)

can.cd()
hAxisBias.Draw()
lineAtZero.Draw("same")
grNomBias.Draw("P")
grMixedBias.SetLineColor(ROOT.kBlue)
grMixedBias.SetMarkerColor(ROOT.kBlue)
grMixedBias.Draw("P")
grSigInjBias.SetMarkerColor(ROOT.kRed)
grSigInjBias.Draw("P")

can.SaveAs(path+"biasVsInjected.pdf")
