#!/usr/bin/env python3

import ROOT
from array import array

ROOT.gROOT.SetBatch(True)

# ================================================================
# Input
# ================================================================

file_path = "/storage4/store/user/zhipeng/Data_For_ReWeight_v2/2022E_No1.root"
tree_name = "JpsiTree"

data_file = ROOT.TFile.Open(file_path)

if not data_file or data_file.IsZombie():
    raise RuntimeError("Cannot open input ROOT file")

input_tree = data_file.Get(tree_name)

if not input_tree:
    raise RuntimeError(f"Cannot find tree: {tree_name}")

# ================================================================
# Variables
# ================================================================

MASS_MIN = 2.5
MASS_MAX = 3.5

CTAU_MIN = -0.03
CTAU_MAX = 0.12

mass = ROOT.RooRealVar("mass","M(J/#psi)",MASS_MIN,MASS_MAX)

ctau = ROOT.RooRealVar( "ctau","c#tau(J/#psi) [cm]",CTAU_MIN,CTAU_MAX)

# ================================================================
# pT and |y| binning
# ================================================================

yBins = array("d",[ 0.0, 0.3, 0.6, 0.9, 1.2, 1.8,2.4 ])

ptBins = array( "d", [10.0,10.5,11.0,11.5, 12.0,13.0, 14.0, 15.0,   16.0,17.0, 18.0, 19.0, 20.0,21.0,22.0,23.0,24.0, 25.0, 26.0,  27.0,28.0,29.0,30.0,32.0,34.0, 36.0, 38.0, 42.0,46.0,50.0,60.0,75.0,95.0,120.0,200.0])

pt = ROOT.RooRealVar("pt","p_{T}(J/#psi) [GeV]", ptBins[0],ptBins[-1])

rapidity = ROOT.RooRealVar(
    "rapidity",
    "|y(J/#psi)|",
    yBins[0],
    yBins[-1]
)

# ================================================================
# RooDataSet
# ================================================================

data = ROOT.RooDataSet(
    "data",
    "data",
    ROOT.RooArgSet(
        mass,
        ctau,
        pt,
        rapidity
    ),
    ROOT.RooFit.Import(input_tree)
)

print("==========================================")
print("Input data entries:", data.numEntries())
print("==========================================")

# ================================================================
# Mass PDF
# Signal = Double Gaussian + Double CB
# ================================================================

mean = ROOT.RooRealVar(
    "mean",
    "m_{J/#psi}",
    3.0964,
    3.08,
    3.11
)

# ------------------------------------------------
# Double Gaussian
# ------------------------------------------------

gSigma1 = ROOT.RooRealVar(
    "gSigma1",
    "gSigma1",
    0.015,
    0.010,
    0.030
)

gSigma2 = ROOT.RooRealVar(
    "gSigma2",
    "gSigma2",
    0.023,
    0.013,
    0.050
)

pdf_gauss1 = ROOT.RooGaussian(
    "pdf_gauss1",
    "Gaussian 1",
    mass,
    mean,
    gSigma1
)

pdf_gauss2 = ROOT.RooGaussian(
    "pdf_gauss2",
    "Gaussian 2",
    mass,
    mean,
    gSigma2
)

fracG1 = ROOT.RooRealVar("fracG1", "f_{G,1}", 0.6,0.0,  1.0)

pdf_doubleGauss = ROOT.RooAddPdf(
    "pdf_doubleGauss",
    "Double Gaussian",
    ROOT.RooArgList(
        pdf_gauss1,
        pdf_gauss2
    ),
    ROOT.RooArgList(
        fracG1
    )
)

# ------------------------------------------------
# Double CB
# ------------------------------------------------

cbSigmaL = ROOT.RooRealVar(
    "cbSigmaL",
    "sigma Left",
    0.016,
    0.001,
    0.035
)

cbAlphaL = ROOT.RooRealVar(
    "cbAlphaL",
    "alpha Left",
    1.5,
    0.2,
    5.0
)

cbNL = ROOT.RooRealVar(
    "cbNL",
    "n Left",
    3.5,
    1.1,
    20.0
)

pdf_cb_left = ROOT.RooCBShape(
    "pdf_cb_left",
    "Left Tail CB",
    mass,
    mean,
    cbSigmaL,
    cbAlphaL,
    cbNL
)

cbSigmaR = ROOT.RooRealVar(
    "cbSigmaR",
    "sigma Right",
    0.016,
    0.01,
    0.050
)

cbAlphaR = ROOT.RooRealVar(
    "cbAlphaR",
    "alpha Right",
    1.5,
    0.2,
    5.0
)

cbNR = ROOT.RooRealVar("cbNR","n Right", 3.5, 1.1,20.0)

pdf_cb_right = ROOT.RooCBShape(
    "pdf_cb_right",
    "Right Tail CB",
    mass,
    mean,
    cbSigmaR,
    cbAlphaR,
    cbNR
)

fracCB1 = ROOT.RooRealVar(
    "fracCB1",
    "fracCB1",
    0.5,
    0.0,
    1.0
)

pdf_doubleCB = ROOT.RooAddPdf(
    "pdf_doubleCB",
    "Double CB",
    ROOT.RooArgList(
        pdf_cb_left,
        pdf_cb_right
    ),
    ROOT.RooArgList(
        fracCB1
    )
)

# ------------------------------------------------
# Total signal mass
# ------------------------------------------------

fracGaussInSig = ROOT.RooRealVar(
    "fracGaussInSig",
    "Gaussian fraction in signal",
    0.8,
    0.0,
    1.0
)

pdf_sigMass = ROOT.RooAddPdf(
    "pdf_sigMass",
    "J/#psi signal mass",
    ROOT.RooArgList(
        pdf_doubleGauss,
        pdf_doubleCB
    ),
    ROOT.RooArgList(
        fracGaussInSig
    )
)

# ================================================================
# Mass background
# ================================================================

chebP1 = ROOT.RooRealVar(
    "chebP1",
    "chebP1",
    -0.66076225,
    -1.0,
    1.0
)

chebP2 = ROOT.RooRealVar(
    "chebP2",
    "chebP2",
    0.030170337,
    -1.0,
    1.0
)

chebP3 = ROOT.RooRealVar(
    "chebP3",
    "chebP3",
    0.053274777,
    -1.0,
    1.0
)

pdf_bkgMass = ROOT.RooChebychev(
    "pdf_bkgMass",
    "Background mass",
    mass,
    ROOT.RooArgList(
        chebP1,
        chebP2,
        chebP3
    )
)

# ================================================================
# Prompt J/psi ctau
# Nested 3-Gaussian
# ================================================================

ctRes_mean = ROOT.RooRealVar(
    "ctRes_mean",
    "ctRes_mean",
    0.0,
    -0.001,
    0.001
)

ctRes_sig1 = ROOT.RooRealVar(
    "ctRes_sig1",
    "ctRes_sig1",
    2.20688e-03,
    0.0005,
    0.050
)

ctRes_sig2 = ROOT.RooRealVar(
    "ctRes_sig2",
    "ctRes_sig2",
    1.53310e-03,
    0.0005,
    0.010
)

ctRes_sig3 = ROOT.RooRealVar(
    "ctRes_sig3",
    "ctRes_sig3",
    1.53310e-03,
    0.0005,
    0.010
)

ctRes_frac1 = ROOT.RooRealVar(
    "ctRes_frac1",
    "ctRes_frac1",
    8.98616e-01,
    0.0,
    1.0
)

ctRes_frac2 = ROOT.RooRealVar(
    "ctRes_frac2",
    "ctRes_frac2",
    8.98616e-01,
    0.0,
    1.0
)

pdf_ctPrompt_g1 = ROOT.RooGaussian(
    "pdf_ctPrompt_g1",
    "Prompt Gaussian 1",
    ctau,
    ctRes_mean,
    ctRes_sig1
)

pdf_ctPrompt_g2 = ROOT.RooGaussian(
    "pdf_ctPrompt_g2",
    "Prompt Gaussian 2",
    ctau,
    ctRes_mean,
    ctRes_sig2
)

pdf_ctPrompt_g3 = ROOT.RooGaussian(
    "pdf_ctPrompt_g3",
    "Prompt Gaussian 3",
    ctau,
    ctRes_mean,
    ctRes_sig3
)

pdf_ctPrompt_g12 = ROOT.RooAddPdf(
    "pdf_ctPrompt_g12",
    "Prompt G1+G2",
    ROOT.RooArgList(
        pdf_ctPrompt_g1,
        pdf_ctPrompt_g2
    ),
    ROOT.RooArgList(
        ctRes_frac1
    )
)

pdf_ctPrompt = ROOT.RooAddPdf(
    "pdf_ctPrompt",
    "Prompt J/#psi c#tau",
    ROOT.RooArgList(
        pdf_ctPrompt_g12,
        pdf_ctPrompt_g3
    ),
    ROOT.RooArgList(
        ctRes_frac2
    )
)

# ================================================================
# Non-prompt J/psi
# Single-sided exponential convolved with Double Gaussian
# ================================================================

ctRes_frac3 = ROOT.RooRealVar(
    "ctRes_frac3",
    "ctRes_frac3",
    8.98616e-01,
    0.0,
    1.0
)

ctNonPrompt_tau = ROOT.RooRealVar(
    "ctNonPrompt_tau",
    "Non-prompt lifetime",
    3.53161e-02,
    0.005,
    0.15
)

resModel_g1 = ROOT.RooGaussModel(
    "resModel_g1",
    "Non-prompt resolution 1",
    ctau,
    ctRes_mean,
    ctRes_sig1
)

resModel_g2 = ROOT.RooGaussModel(
    "resModel_g2",
    "Non-prompt resolution 2",
    ctau,
    ctRes_mean,
    ctRes_sig2
)

resModel_sum = ROOT.RooAddModel(
    "resModel_sum",
    "Double-Gaussian resolution",
    ROOT.RooArgList(
        resModel_g1,
        resModel_g2
    ),
    ROOT.RooArgList(
        ctRes_frac3
    )
)

pdf_ctNonPrompt = ROOT.RooDecay(
    "pdf_ctNonPrompt",
    "Non-prompt J/#psi c#tau",
    ctau,
    ctNonPrompt_tau,
    resModel_sum,
    ROOT.RooDecay.SingleSided
)

# ================================================================
# Background ctau
# ================================================================

ctBkg_mean = ROOT.RooRealVar(
    "ctBkg_mean",
    "Background c#tau mean",
    0.0,
    -0.002,
    0.002
)

ctBkg_sig1 = ROOT.RooRealVar(
    "ctBkg_sig1",
    "Background c#tau sigma 1",
    0.0015,
    0.0005,
    0.0035
)

ctBkg_sig2 = ROOT.RooRealVar(
    "ctBkg_sig2",
    "Background c#tau sigma 2",
    0.0050,
    0.0036,
    0.0200
)

ctBkg_frac1 = ROOT.RooRealVar(
    "ctBkg_frac1",
    "Background Gaussian fraction",
    0.75,
    0.0,
    1.0
)

resModelBkg_g1 = ROOT.RooGaussModel(
    "resModelBkg_g1",
    "Bkg Gaussian resolution 1",
    ctau,
    ctBkg_mean,
    ctBkg_sig1
)

resModelBkg_g2 = ROOT.RooGaussModel(
    "resModelBkg_g2",
    "Bkg Gaussian resolution 2",
    ctau,
    ctBkg_mean,
    ctBkg_sig2
)

resModelBkg_sum = ROOT.RooAddModel(
    "resModelBkg_sum",
    "Bkg Double-Gaussian resolution",
    ROOT.RooArgList(
        resModelBkg_g1,
        resModelBkg_g2
    ),
    ROOT.RooArgList(
        ctBkg_frac1
    )
)

pdf_ctBkg_prompt = ROOT.RooAddModel(
    "pdf_ctBkg_prompt",
    "Bkg prompt component",
    ROOT.RooArgList(
        resModelBkg_g1,
        resModelBkg_g2
    ),
    ROOT.RooArgList(
        ctBkg_frac1
    )
)

ctBkg_tau_right = ROOT.RooRealVar(
    "ctBkg_tau_right",
    "Bkg right lifetime",
    0.03,
    0.005,
    0.2
)

pdf_ctBkg_tailRight = ROOT.RooDecay(
    "pdf_ctBkg_tailRight",
    "Bkg right tail",
    ctau,
    ctBkg_tau_right,
    resModelBkg_sum,
    ROOT.RooDecay.SingleSided
)

ctBkg_tau_left = ROOT.RooRealVar(
    "ctBkg_tau_left",
    "Bkg left lifetime",
    0.005,
    0.001,
    0.02
)

pdf_ctBkg_tailLeft = ROOT.RooDecay(
    "pdf_ctBkg_tailLeft",
    "Bkg left tail",
    ctau,
    ctBkg_tau_left,
    resModelBkg_sum,
    ROOT.RooDecay.Flipped
)

ctBkg_f_prompt = ROOT.RooRealVar(
    "ctBkg_f_prompt",
    "Bkg prompt fraction",
    0.4,
    0.0,
    1.0
)

ctBkg_f_right = ROOT.RooRealVar(
    "ctBkg_f_right",
    "Bkg right-tail fraction",
    0.5,
    0.0,
    1.0
)

pdf_ctBkg = ROOT.RooAddPdf(
    "pdf_ctBkg",
    "Final background c#tau PDF",
    ROOT.RooArgList(
        pdf_ctBkg_prompt,
        pdf_ctBkg_tailRight,
        pdf_ctBkg_tailLeft
    ),
    ROOT.RooArgList(
        ctBkg_f_prompt,
        ctBkg_f_right
    )
)

# ================================================================
# Final 2D PDF
# ================================================================

pdf_promptJpsi = ROOT.RooProdPdf(
    "pdf_promptJpsi",
    "Prompt J/#psi (mass x ctau)",
    ROOT.RooArgList(
        pdf_sigMass,
        pdf_ctPrompt
    )
)

pdf_nonPromptJpsi = ROOT.RooProdPdf(
    "pdf_nonPromptJpsi",
    "Non-prompt J/#psi (mass x ctau)",
    ROOT.RooArgList(
        pdf_sigMass,
        pdf_ctNonPrompt
    )
)

pdf_bkg = ROOT.RooProdPdf(
    "pdf_bkg",
    "Combinatorial background (mass x ctau)",
    ROOT.RooArgList(
        pdf_bkgMass,
        pdf_ctBkg
    )
)

# ================================================================
# Extended yields
# ================================================================

totalEntries = data.numEntries()

nPrompt = ROOT.RooRealVar(
    "nPrompt",
    "prompt yield",
    totalEntries * 0.5,
    0.0,
    totalEntries
)

nNonPrompt = ROOT.RooRealVar(
    "nNonPrompt",
    "non-prompt yield",
    totalEntries * 0.3,
    0.0,
    totalEntries
)

nBkg = ROOT.RooRealVar(
    "nBkg",
    "background yield",
    totalEntries * 0.2,
    0.0,
    totalEntries
)

pdf_total = ROOT.RooAddPdf(
    "pdf_total",
    "Total 2D PDF",
    ROOT.RooArgList(
        pdf_promptJpsi,
        pdf_nonPromptJpsi,
        pdf_bkg
    ),
    ROOT.RooArgList(
        nPrompt,
        nNonPrompt,
        nBkg
    )
)

# ================================================================
# Your latest fit parameters
# ================================================================

cbAlphaL.setVal(4.9770999)
cbAlphaR.setVal(0.32946295)

cbNL.setVal(1.1001213)
cbNR.setVal(1.1001162)

cbSigmaL.setVal(0.0050181123)
cbSigmaR.setVal(0.040315067)

chebP1.setVal(-0.1210705)
chebP2.setVal(-0.098898918)
chebP3.setVal(-0.27450534)

ctBkg_f_prompt.setVal(0.1497091)
ctBkg_f_right.setVal(0.73010626)

ctBkg_frac1.setVal(0.39952204)

ctBkg_mean.setVal(6.5587646e-05)

ctBkg_sig1.setVal(0.0018069213)
ctBkg_sig2.setVal(0.0036045979)

ctBkg_tau_left.setVal(0.012588601)
ctBkg_tau_right.setVal(0.065129147)

ctNonPrompt_tau.setVal(0.037107806)

ctRes_frac1.setVal(0.60460739)
ctRes_frac2.setVal(0.9411862)
ctRes_frac3.setVal(0.34224434)

ctRes_mean.setVal(4.5059156e-05)

ctRes_sig1.setVal(0.00080095692)
ctRes_sig2.setVal(0.001773384)
ctRes_sig3.setVal(0.00052525272)

fracCB1.setVal(0.034625593)
fracG1.setVal(0.68873064)

fracGaussInSig.setVal(0.78135502)

gSigma1.setVal(0.039894428)
gSigma2.setVal(0.025570624)

mean.setVal(3.0945947)

# ================================================================
# Initial yields
# ================================================================

nBkg.setVal(21298.862)
nNonPrompt.setVal(31862.065)
nPrompt.setVal(46838.984)

# ================================================================
# Fit
# ================================================================

print("==========================================")
print("Starting 2D mass + ctau fit")
print("==========================================")

fitResult = pdf_total.fitTo(
    data,
    ROOT.RooFit.Extended(ROOT.kTRUE),
    ROOT.RooFit.Save(ROOT.kTRUE),
    ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
    ROOT.RooFit.Strategy(0),
    ROOT.RooFit.NumCPU(4)
)

print("==========================================")
print("Fit result")
print("==========================================")

print("EDM:", fitResult.edm())
print("Fit status:", fitResult.status())
print("Covariance quality:", fitResult.covQual())

fitResult.Print()

print("==========================================")
print("Final parameters")
print("==========================================")

for par in fitResult.floatParsFinal():
    print(
        f"{par.GetName()}.setVal({par.getVal():.8g})"
    )

# ================================================================
# sPlot
# ================================================================

print("==========================================")
print("Running sPlot")
print("==========================================")

sData = ROOT.RooStats.SPlot(
    "sData",
    "sData",
    data,
    pdf_total,
    ROOT.RooArgList(
        nPrompt,
        nNonPrompt,
        nBkg
    )
)

# ================================================================
# Prompt J/psi pT - |y| distribution
# ================================================================

hPromptPtY = ROOT.TH2D(
    "hPromptPtY",
    "Prompt J/#psi;p_{T}(J/#psi) [GeV];|y(J/#psi)|",
    len(ptBins) - 1,
    ptBins,
    len(yBins) - 1,
    yBins
)

for i in range(data.numEntries()):

    row = data.get(i)

    thisPt = row.getRealValue("pt")
    thisY = abs(row.getRealValue("rapidity"))

    promptWeight = row.getRealValue("nPrompt_sw")

    if (
        thisPt < ptBins[0]
        or thisPt >= ptBins[-1]
    ):
        continue

    if (
        thisY < yBins[0]
        or thisY >= yBins[-1]
    ):
        continue

    hPromptPtY.Fill(
        thisPt,
        thisY,
        promptWeight
    )

# ================================================================
# Projection
# ================================================================

hPromptPt = hPromptPtY.ProjectionX(
    "hPromptPt"
)

hPromptY = hPromptPtY.ProjectionY(
    "hPromptY"
)

print("==========================================")
print("Prompt sPlot yield")
print("==========================================")

print(
    "2D prompt yield =",
    hPromptPtY.Integral()
)

print(
    "pT projection =",
    hPromptPt.Integral()
)

print(
    "|y| projection =",
    hPromptY.Integral()
)

print(
    "Fitted nPrompt =",
    nPrompt.getVal()
)

# ================================================================
# Fit cross-check: Mass
# ================================================================

plotMass = mass.frame(
    ROOT.RooFit.Name("plotMass"),
    ROOT.RooFit.Title("J/#psi mass fit")
)

data.plotOn(
    plotMass,
    ROOT.RooFit.Binning(100)
)

pdf_total.plotOn(
    plotMass,
    ROOT.RooFit.LineColor(ROOT.kRed),
    ROOT.RooFit.LineWidth(2)
)

pdf_total.plotOn(
    plotMass,
    ROOT.RooFit.Components("pdf_promptJpsi"),
    ROOT.RooFit.LineColor(ROOT.kGreen + 2),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

pdf_total.plotOn(
    plotMass,
    ROOT.RooFit.Components("pdf_nonPromptJpsi"),
    ROOT.RooFit.LineColor(ROOT.kBlue),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

pdf_total.plotOn(
    plotMass,
    ROOT.RooFit.Components("pdf_bkg"),
    ROOT.RooFit.LineColor(ROOT.kMagenta),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

# ================================================================
# Fit cross-check: Ctau
# ================================================================

plotCtau = ctau.frame(
    ROOT.RooFit.Name("plotCtau"),
    ROOT.RooFit.Title("J/#psi c#tau fit")
)

data.plotOn(
    plotCtau,
    ROOT.RooFit.Binning(100)
)

pdf_total.plotOn(
    plotCtau,
    ROOT.RooFit.LineColor(ROOT.kRed),
    ROOT.RooFit.LineWidth(2)
)

pdf_total.plotOn(
    plotCtau,
    ROOT.RooFit.Components("pdf_promptJpsi"),
    ROOT.RooFit.LineColor(ROOT.kGreen + 2),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

pdf_total.plotOn(
    plotCtau,
    ROOT.RooFit.Components("pdf_nonPromptJpsi"),
    ROOT.RooFit.LineColor(ROOT.kBlue),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

pdf_total.plotOn(
    plotCtau,
    ROOT.RooFit.Components("pdf_bkg"),
    ROOT.RooFit.LineColor(ROOT.kMagenta),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

# ================================================================
# Save
# ================================================================

output_file = ROOT.TFile(
    "fitOutput.root",
    "RECREATE"
)

hPromptPtY.Write()
hPromptPt.Write()
hPromptY.Write()

plotMass.Write()
plotCtau.Write()

fitResult.Write("fitResult")

output_file.Close()

data_file.Close()

print("==========================================")
print("Saved:")
print("  fitOutput.root")
print("  hPromptPtY")
print("  hPromptPt")
print("  hPromptY")
print("==========================================")
