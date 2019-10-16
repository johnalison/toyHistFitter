import ROOT
rand = ROOT.TRandom3()
ROOT.gROOT.SetBatch(True)


hBkg = ROOT.TH1F("Background","backgound",100,-10,10);
hSig = ROOT.TH1F("Signal","signal",100,-10,10);
hSig.SetLineColor(ROOT.kRed)

for i in range(10000):
    hBkg.Fill(rand.Gaus(0,1))
    hSig.Fill(rand.Gaus(4,1))


def makeFakeData(bkg,sig,mu):
    hFakeData = bkg.Clone()

    hFakeData.Add(sig, mu)
    hFakeData.SetName("FakeData")
    hFakeData.SetTitle("FakeData")
    
    return hFakeData


def sPlusB_func(x,par):    
    nBkg = hBkg.GetBinContent(hBkg.FindBin(x[0]))
    nSig = hSig.GetBinContent(hSig.FindBin(x[0]))
    return nBkg+par[0]*nSig


tf1_sPlusB_func = ROOT.TF1("tf1_sPlusB",sPlusB_func,-10,10,1)
tf1_sPlusB_func.SetParameter(0,1)
print tf1_sPlusB_func.Eval(5)



can = ROOT.TCanvas("canSB","canSB")
can.cd()
hBkg.Draw()
hSig.Draw("same")
can.SaveAs("RawSandB.pdf")


hFake = makeFakeData(hBkg,hSig,0.5)
hFake.Draw()
can.SaveAs("FakeData.pdf")


hFake.Fit(tf1_sPlusB_func,"Q")
print "fitted Value",tf1_sPlusB_func.GetParameter(0)
can.SaveAs("FakeDataWithFit.pdf")

hFitAxis = ROOT.TH1F("axis","Injected vs Fitted mu;injected mu;fitted mu",1,0,1)
hFitAxis.SetStats(0)
fitResults = ROOT.TGraph(10)
fitResults.SetMarkerStyle(20)
fitResults.SetMarkerSize(1.2)
iPt = 0

for i in range(0,100,10):
    mu = float(i)/100
    hFake = makeFakeData(hBkg,hSig,mu)
    tf1_sPlusB_func.SetParameter(0,1)
    hFake.Fit(tf1_sPlusB_func,"Q")
    mu_fit = tf1_sPlusB_func.GetParameter(0)
    iPt += 1
    print "True mu",mu,"Fitted mu",mu_fit
    fitResults.SetPoint(iPt,mu,mu_fit)

can.cd()
hFitAxis.Draw()
fitResults.Draw("PL")
can.SaveAs("FitResults.pdf")
