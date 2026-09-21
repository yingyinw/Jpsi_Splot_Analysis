#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ROOT
import glob
from array import array
from DataFormats.FWLite import Events, Handle

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# ============================================================
# Input
# ============================================================

refFile = "fitOutput_final.root"

files = glob.glob("/eos/user/z/zhipeng/MC_Particle_Gun_For_Acc/results/step0_*.root")

if len(files) == 0:
    raise RuntimeError("No input Particle Gun ROOT files found")

print("==============================")
print("Number of input files =", len(files))
print("==============================")

# ============================================================
# Load SPlot prompt reference
# ============================================================

fRef = ROOT.TFile.Open(refFile)

if not fRef or fRef.IsZombie():
    raise RuntimeError("Cannot open fitOutput_final.root")

hRef = fRef.Get("hPromptPtY")

if not hRef:
    raise RuntimeError("Cannot find hPromptPtY in fitOutput_final.root")

print("==============================")
print("Reference histogram")
print("Name      =", hRef.GetName())
print("Integral  =", hRef.Integral())
print("Nbins X   =", hRef.GetNbinsX())
print("Nbins Y   =", hRef.GetNbinsY())
print("==============================")

# ============================================================
# GEN particle handle
# ============================================================

handle = Handle("std::vector<reco::GenParticle>")
label = ("genParticles",)

# ============================================================
# Analysis range 
# ============================================================ 

xMin = hRef.GetXaxis().GetXmin() 
xMax = hRef.GetXaxis().GetXmax() 

yMin = hRef.GetYaxis().GetXmin() 
yMax = hRef.GetYaxis().GetXmax()

# ============================================================
# First pass
# Build MC generated GEN J/psi pT-|y| distribution
# ============================================================

hGenPtY = hRef.Clone("hGenPtY")
hGenPtY.Reset()
hGenPtY.SetTitle("MC generated J/#psi;p_{T} [GeV];|y|")

events = Events(files)

nEvents = 0
nJpsi = 0

for iev, event in enumerate(events):

    nEvents += 1

    if iev % 10000 == 0:
        print("First pass:", iev)

    event.getByLabel(label, handle)
    genParticles = handle.product()

    for p in genParticles:

        if p.pdgId() != 443:
            continue

        if not p.isLastCopy():
            continue

        nJpsi += 1

        pt = p.pt()
        y = abs(p.rapidity())

        if not (xMin <= pt < xMax):
            continue

        if not (yMin <= y < yMax):
            continue

        hGenPtY.Fill(pt, y)

print("==============================")
print("First pass finished")
print("Events      =", nEvents)
print("J/psi       =", nJpsi)
print("MC integral =", hGenPtY.Integral())
print("==============================")

# ============================================================
# Build pT-|y| weight map
#
# weight =
#     normalized prompt DATA
#     ----------------------
#     normalized Particle Gun MC
# ============================================================

hWeight = hRef.Clone("hPtYWeight")
hWeight.Reset()

hWeight.SetTitle("p_{T}-|y| reweighting factor;p_{T} [GeV];|y|")

refIntegral = hRef.Integral()
genIntegral = hGenPtY.Integral()

if refIntegral <= 0:
    raise RuntimeError(
        "Reference histogram hPromptPtY has zero integral")

if genIntegral <= 0:
    raise RuntimeError(
        "MC histogram hGenPtY has zero integral")

for ix in range(1, hRef.GetNbinsX() + 1):

    for iy in range(1, hRef.GetNbinsY() + 1):

        refContent = hRef.GetBinContent(ix, iy)
        genContent = hGenPtY.GetBinContent(ix, iy)

        refFrac = refContent / refIntegral
        genFrac = genContent / genIntegral

        if genFrac > 0:
            weight = refFrac / genFrac
        else:
            weight = 0.0

        hWeight.SetBinContent(ix, iy, weight)

print("==============================")
print("Weight map created")
print("Weight minimum =", hWeight.GetMinimum())
print("Weight maximum =", hWeight.GetMaximum())
print("==============================")

# ============================================================
# Weight function
# ============================================================

def getPtYWeight(pt, y):

    y = abs(y)

    if not (xMin <= pt < xMax):
        return 0.0

    if not (yMin <= y < yMax):
        return 0.0

    ix = hWeight.GetXaxis().FindFixBin(pt)
    iy = hWeight.GetYaxis().FindFixBin(y)

    if ix < 1 or ix > hWeight.GetNbinsX():
        return 0.0

    if iy < 1 or iy > hWeight.GetNbinsY():
        return 0.0

    weight = hWeight.GetBinContent(ix, iy)

    if weight <= 0:
        return 0.0

    return weight

# ============================================================
# Unweighted acceptance histograms
# ============================================================

hDenomOriginal = hRef.Clone("hDenomOriginal")
hDenomOriginal.Reset()

hDenomOriginal.SetTitle("Unweighted GEN J/#psi denominator;""p_{T} [GeV];|y|")

hNumerOriginal = hRef.Clone("hNumerOriginal")
hNumerOriginal.Reset()

hNumerOriginal.SetTitle("Unweighted accepted J/#psi numerator;""p_{T} [GeV];|y|")

# ============================================================
# Weighted acceptance Histograms
# ============================================================

hDenomWeighted = hRef.Clone("hDenomWeighted")
hDenomWeighted.Reset()

hDenomWeighted.SetTitle("Weighted GEN J/#psi denominator;" "p_{T} [GeV];|y|")

hNumerWeighted = hRef.Clone("hNumerWeighted")
hNumerWeighted.Reset()

hNumerWeighted.SetTitle("Weighted accepted J/#psi numerator;" "p_{T} [GeV];|y|")

hPtOriginal = ROOT.TH1D("hPtOriginal","GEN J/#psi p_{T};p_{T} [GeV];Events",40,xMin,xMax)

hPtWeighted = ROOT.TH1D("hPtWeighted","Weighted GEN J/#psi p_{T};p_{T} [GeV];Weighted events", 40,xMin,xMax)

hYOriginal = ROOT.TH1D("hYOriginal","GEN J/#psi |y|;|y|;Events",30,yMin, yMax)

hYWeighted = ROOT.TH1D("hYWeighted", "Weighted GEN J/#psi |y|;|y|;Weighted events",30,yMin,yMax)

for h in [
    hDenomOriginal,
    hNumerOriginal,
    hDenomWeighted,
    hNumerWeighted,
    hPtOriginal,
    hPtWeighted,
    hYOriginal,
    hYWeighted
]:
    h.Sumw2()

# ============================================================
# Muon acceptance
# ============================================================

def muon_pass_criteria(mu):

    if abs(mu.eta()) > 2.4:
        return False

    if mu.pt() < 5.0:
        return False

    return True

# ============================================================
# Find final-state muons
# ============================================================

def find_final_muons(p):

    muons = []

    if abs(p.pdgId()) == 13:

        if p.numberOfDaughters() == 0:

            muons.append(p)

        else:

            for i in range(p.numberOfDaughters()):

                muons.extend(
                    find_final_muons(
                        p.daughter(i)
                    )
                )

        return muons

    for i in range(p.numberOfDaughters()):
        muons.extend(
            find_final_muons(p.daughter(i))
        )

    return muons

# ============================================================
# Second pass
# Apply pT-|y| weight and calculate acceptance
# ============================================================

events = Events(files)

nSelected = 0
nAccepted = 0
sumWeights = 0.0

for iev, event in enumerate(events):

    if iev % 10000 == 0:
        print("Second pass:", iev)

    event.getByLabel(label, handle)

    genParticles = handle.product()

    for p in genParticles:

        if p.pdgId() != 443:
            continue
            
        if not p.isLastCopy():
            continue

        pt = p.pt()
        y = abs(p.rapidity())

        if not (xMin <= pt < xMax):
            continue

        if not (yMin <= y < yMax):
            continue

        weight = getPtYWeight(pt, y)

        if weight <= 0:
            continue

        nSelected += 1

        sumWeights += weight

        hPtOriginal.Fill(pt)
        hPtWeighted.Fill(pt, weight)

        hYOriginal.Fill(y)
        hYWeighted.Fill(y, weight)
       
        hDenomOriginal.Fill(pt, y)
        hDenomWeighted.Fill(pt,y,weight)
      
        muons = find_final_muons(p)

        if len(muons) != 2:
            continue

        if (
            muon_pass_criteria(muons[0])
            and
            muon_pass_criteria(muons[1])
        ):
            hNumerOriginal.Fill(pt, y)
            hNumerWeighted.Fill(pt, y,weight)

            nAccepted += 1

# ============================================================
# Unweighted acceptance
# ============================================================

hAcceptanceOriginal = hNumerOriginal.Clone("hAcceptanceOriginal")

hAcceptanceOriginal.SetTitle("Unweighted J/#psi acceptance;""p_{T} [GeV];|y|")

hAcceptanceOriginal.Divide(hNumerOriginal,hDenomOriginal,1.0,1.0,"")

# ============================================================
# Weighted acceptance
# ============================================================

hAcceptanceWeighted = hNumerWeighted.Clone("hAcceptanceWeighted")

hAcceptanceWeighted.SetTitle( "Weighted J/#psi acceptance;" "p_{T} [GeV];|y|")

hAcceptanceWeighted.Divide(hNumerWeighted,hDenomWeighted,1.0,1.0,"")

# ============================================================
# Final checks
# ============================================================

denomIntegral = hDenomWeighted.Integral()
numerIntegral = hNumerWeighted.Integral()

# ------------------------------------------------------------
# Calculate total statistical errors
# For weighted histograms:
# sigma_total^2 = sum_i (sigma_i)^2
# ------------------------------------------------------------

denomError2 = 0.0
numerError2 = 0.0

for ix in range(1, hDenomWeighted.GetNbinsX() + 1):

    for iy in range(1, hDenomWeighted.GetNbinsY() + 1):

        denomBinError = hDenomWeighted.GetBinError(ix, iy)
        numerBinError = hNumerWeighted.GetBinError(ix, iy)

        denomError2 += denomBinError * denomBinError
        numerError2 += numerBinError * numerBinError

denomError = denomError2 ** 0.5
numerError = numerError2 ** 0.5

# ------------------------------------------------------------
# Integrated acceptance
# A = N / D
# ------------------------------------------------------------

if denomIntegral > 0 and numerIntegral > 0:

    integratedAcceptance = (numerIntegral / denomIntegral)

    # --------------------------------------------------------
    # Error propagation
    #
    # sigma_A = A * sqrt( (sigma_N / N)^2+  (sigma_D / D)^2)
    # --------------------------------------------------------

    integratedAcceptanceError = (
        integratedAcceptance*((numerError / numerIntegral) ** 2+ (denomError / denomIntegral) ** 2) ** 0.5 )

else:

    integratedAcceptance = 0.0
    integratedAcceptanceError = 0.0

# ------------------------------------------------------------
# Print final results
# ------------------------------------------------------------

print("==============================")
print("Denominator = {:.8f} +/- {:.8f}".format(denomIntegral,denomError))
print("Numerator   = {:.8f} +/- {:.8f}".format(numerIntegral, numerError))
print("Acceptance  = {:.8f} +/- {:.8f}".format(integratedAcceptance,integratedAcceptanceError))
print("Sum weights = {:.8f}".format(sumWeights))
print("Selected    = {}".format( nSelected))
print("Accepted    = {}".format( nAccepted))
print("==============================")

# ============================================================
# Write acceptance results to TXT
# ============================================================

txtFile = open("Acceptance_reweighted.txt","w")
#txtFile = open("Acceptance_reweighted_test.txt","w")

txtFile.write("============================================================\n")
txtFile.write("Weighted J/psi Acceptance\n")
txtFile.write("============================================================\n")
txtFile.write("Number of input files       = {}\n".format(len(files)))
txtFile.write("Total events                 = {}\n".format(nEvents))
txtFile.write("Total GEN J/psi              = {}\n".format(nJpsi))
txtFile.write("Selected J/psi               = {}\n".format(nSelected))
txtFile.write("Accepted J/psi               = {}\n".format(nAccepted ))
txtFile.write("Sum of weights               = {:.8f}\n".format(sumWeights))
txtFile.write("Integrated acceptance        = {:.8f} +/- {:.8f}\n".format(integratedAcceptance, integratedAcceptanceError))
txtFile.write("\n")

# ============================================================
# Bin-by-bin unweighted and weighted acceptance
# ============================================================

txtFile.write("============================================================\n")
txtFile.write("Bin-by-bin acceptance comparison\n")
txtFile.write("============================================================\n")

txtFile.write(
    "{:<10} {:<10} {:<10} {:<10} "
    "{:<18} {:<18} {:<18} "
    "{:<18} {:<18} {:<18}\n".format(
        "pT_low",
        "pT_high",
        "y_low",
        "y_high",
        "Denominator",
        "Numerator",
        "Acceptance",
        "Weighted_Denominator",
        "Weighted_Numerator",
        "Weighted_Acceptance"
    )
)

txtFile.write("------------------------------------------------------------------------------------------------------\n")

for ix in range(
    1,
    hAcceptanceWeighted.GetNbinsX() + 1
):

    ptLow = hAcceptanceWeighted.GetXaxis().GetBinLowEdge(ix)
    ptHigh = hAcceptanceWeighted.GetXaxis().GetBinUpEdge(ix)

    for iy in range(
        1,
        hAcceptanceWeighted.GetNbinsY() + 1
    ):

        yLow = hAcceptanceWeighted.GetYaxis().GetBinLowEdge(iy)
        yHigh = hAcceptanceWeighted.GetYaxis().GetBinUpEdge(iy)

        # ----------------------------------------------------
        # Unweighted
        # ----------------------------------------------------

        denominatorOriginal = (hDenomOriginal.GetBinContent(ix, iy))

        numeratorOriginal = (hNumerOriginal.GetBinContent(ix, iy))

        acceptanceOriginal = (hAcceptanceOriginal.GetBinContent(ix, iy))

        # ----------------------------------------------------
        # Weighted
        # ----------------------------------------------------

        denominatorWeighted = (hDenomWeighted.GetBinContent(ix, iy))
        
        numeratorWeighted = (hNumerWeighted.GetBinContent(ix, iy))

        acceptanceWeighted = (hAcceptanceWeighted.GetBinContent(ix, iy))

        txtFile.write(
            "{:<10.3f} {:<10.3f} "
            "{:<10.3f} {:<10.3f} "
            "{:<18.8f} {:<18.8f} {:<18.8f} "
            "{:<18.8f} {:<18.8f} {:<18.8f}\n".format(
                ptLow,
                ptHigh,
                yLow,
                yHigh,
                denominatorOriginal,
                numeratorOriginal,
                acceptanceOriginal,
                denominatorWeighted,
                numeratorWeighted,
                acceptanceWeighted
            )
        )

txtFile.write(
    "============================================================\n")

print("==============================")
print("TXT output:")
print("Acceptance_reweighted.txt")
#print("Acceptance_reweighted_test.txt")
print("==============================")

# ============================================================
# Draw validation plots
# ============================================================

ROOT.gStyle.SetOptStat(1110)

def draw_compare(h1, h2, name):

    canvas = ROOT.TCanvas(name,name,900,700)

    h1.SetLineColor(ROOT.kBlack)
    h1.SetLineWidth(2)

    h2.SetLineColor(ROOT.kRed)
    h2.SetLineWidth(2)

    ymax = max(h1.GetMaximum(),h2.GetMaximum())

    h1.SetMaximum(ymax * 1.35)

    h1.Draw("HIST")
    h2.Draw("HIST SAME")

    canvas.Update()

    # --------------------------------------------------------
    # Move statistics boxes
    # --------------------------------------------------------

    stats1 = h1.FindObject("stats")

    if stats1:

        stats1.SetX1NDC(0.68)
        stats1.SetX2NDC(0.90)
        stats1.SetY1NDC(0.78)
        stats1.SetY2NDC(0.92)

    stats2 = h2.FindObject("stats")

    if stats2:

        stats2.SetX1NDC(0.68)
        stats2.SetX2NDC(0.90)
        stats2.SetY1NDC(0.60)
        stats2.SetY2NDC(0.74)

    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

    leg = ROOT.TLegend(0.18,0.75,0.48, 0.88 )
    leg.SetBorderSize(0)
    leg.AddEntry(h1,"Particle gun","l")
    leg.AddEntry(h2,"After p_{T}-|y| weighting","l")
    leg.Draw()
    canvas.Update()
    canvas.SaveAs(name + ".png" )

# ============================================================
# pT comparison
# ============================================================

draw_compare(hPtOriginal,hPtWeighted,"JpsiPt_before_after")
#draw_compare(hPtOriginal,hPtWeighted,"JpsiPt_before_after_test")

# ============================================================
# |y| comparison
# ============================================================

draw_compare(hYOriginal,hYWeighted,"JpsiY_before_after")
#draw_compare(hYOriginal,hYWeighted,"JpsiY_before_after_test")

# ============================================================
# Draw weight map
# ============================================================

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat("4.3f")

cWeight = ROOT.TCanvas("cWeight","pT-|y| Weight",1000,800)

hWeight.SetTitle("p_{T}-|y| reweighting factor;""p_{T} [GeV];|y|")

hWeight.SetMarkerSize(0.7)
hWeight.SetMinimum(0.0)

hWeight.Draw("COLZ TEXT")
cWeight.Update()

cWeight.SaveAs("Jpsi_pT_y_weight.png")
#cWeight.SaveAs("Jpsi_pT_y_weight_test.png")

# ============================================================
# Create transposed weighted acceptance map
# Original:X = pT Y = |y|
# New:X = |y|  Y = pT
# ============================================================
yBins = array( 'd',[0.0,0.3,0.6,0.9,1.2,1.8, 2.4])

ptBins = array('d',[10,10.5,11,11.5,12,13,14,15, 16,17,18,19, 20,21,22,23,24,25,26,27,28,29,30,32,34,36,38,42,46,50,60,75,95,120,200])

hAcceptanceWeighted_Transpose = ROOT.TH2D(
    "hAcceptanceWeighted_Transpose",
    "Weighted J/#psi acceptance;|y|;p_{T} [GeV]",
    len(yBins) - 1,
    yBins,
    len(ptBins) - 1,
    ptBins)

hAcceptanceWeighted_Transpose.SetDirectory(0)

for ix in range(
    1,
    hAcceptanceWeighted.GetNbinsX() + 1
):

    for iy in range(
        1,
        hAcceptanceWeighted.GetNbinsY() + 1
    ):

        acc = hAcceptanceWeighted.GetBinContent( ix,iy)

        err = hAcceptanceWeighted.GetBinError(ix,iy)

        hAcceptanceWeighted_Transpose.SetBinContent(iy,ix, acc )

        hAcceptanceWeighted_Transpose.SetBinError(iy,ix,err)

# ============================================================
# Draw weighted acceptance map
# X = |y|  Y = pT
# ============================================================

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat("4.2f")

cAcceptanceWeighted = ROOT.TCanvas("cAcceptanceWeighted","Weighted Acceptance Map", 1000,800)

# Logarithmic pT axis
cAcceptanceWeighted.SetLogy()

hAcceptanceWeighted_Transpose.SetTitle("Weighted J/#psi acceptance;""|y|;p_{T} [GeV]")
hAcceptanceWeighted_Transpose.SetMarkerSize(0.7)
hAcceptanceWeighted_Transpose.SetMinimum(0.0)
hAcceptanceWeighted_Transpose.SetMaximum( 1.0)
hAcceptanceWeighted_Transpose.GetXaxis().SetTitle("|y|")
hAcceptanceWeighted_Transpose.GetYaxis().SetTitle("p_{T} [GeV]")
hAcceptanceWeighted_Transpose.GetXaxis().SetRangeUser(0.0,2.4)
hAcceptanceWeighted_Transpose.GetYaxis().SetRangeUser(10.0,200.0)

#hAcceptanceWeighted_Transpose.Draw("COLZ TEXT")
hAcceptanceWeighted_Transpose.Draw("COLZ")
cAcceptanceWeighted.Update()

#cAcceptanceWeighted.SaveAs("Jpsi_acceptance_weighted.png")
cAcceptanceWeighted.SaveAs("Jpsi_acceptance_weighted_v2.png")
#cAcceptanceWeighted.SaveAs("Jpsi_acceptance_weighted_test.png")

# ============================================================
# Save output
# ============================================================

outFile = ROOT.TFile("Acceptance_reweighted.root","RECREATE")
#outFile = ROOT.TFile("Acceptance_reweighted_test.root","RECREATE")

hRef.Write()
hGenPtY.Write()
hWeight.Write()

hDenomWeighted.Write()
hNumerWeighted.Write()

# Original acceptance map
hAcceptanceWeighted.Write()

# Transposed acceptance map
hAcceptanceWeighted_Transpose.Write()

hPtOriginal.Write()
hPtWeighted.Write()

hYOriginal.Write()
hYWeighted.Write()

outFile.Close()
fRef.Close()

# ============================================================
# Finished
# ============================================================

print("==============================")
print("Finished")
print("Events      =", nEvents)
print("J/psi       =", nJpsi)
print("Selected    =", nSelected)
print("Accepted    =", nAccepted)
print("Acceptance  =",integratedAcceptance)
print("==============================")

