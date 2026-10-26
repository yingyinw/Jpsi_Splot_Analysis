
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ROOT
import glob
import os
import math
from array import array
from DataFormats.FWLite import Events, Handle

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(1110)

# ============================================================
# 1. Input
# ============================================================

refFile = "fitOutput_final.root"

inputPattern = (
    "/eos/user/z/zhipeng/MC_Particle_Gun_For_Acc/results/step0_*.root"
)

files = sorted(glob.glob(inputPattern))

if len(files) == 0:
    raise RuntimeError(
        "No input Particle Gun ROOT files found: " + inputPattern
    )

outputRootFileName = "Acceptance_reweighted_v3.root"
outputTextFileName = "Acceptance_reweighted_v3.txt"

print("=" * 70)
print("Input files")
print("=" * 70)
print("Number of files =", len(files))
print("Reference file  =", refFile)
print("=" * 70)

# ============================================================
# 2. Read reference histogram
# ============================================================

fRef = ROOT.TFile.Open(refFile)

if not fRef or fRef.IsZombie():
    raise RuntimeError("Cannot open reference file: " + refFile)

hRef = fRef.Get("hPromptPtY")

if not hRef:
    fRef.Close()
    raise RuntimeError("Cannot find hPromptPtY in " + refFile)

# Detach the histogram from the input file.
hRef.SetDirectory(0)

print("\nReference histogram")
print("Name       =", hRef.GetName())
print("Integral   =", hRef.Integral())
print("Nbins X    =", hRef.GetNbinsX())
print("Nbins Y    =", hRef.GetNbinsY())
print("X range    =", hRef.GetXaxis().GetXmin(),"to",hRef.GetXaxis().GetXmax())
print("Y range    =", hRef.GetYaxis().GetXmin(),"to",hRef.GetYaxis().GetXmax())

# ============================================================
# 3. Fine binning for the weight map
#
# 68 pT bins x 12 |y| bins
# ============================================================

fineYBins = array("d", [
    0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9,
    1.05, 1.2, 1.5, 1.8, 2.1, 2.4
])

finePtBins = array("d", [
    10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
    11.5, 11.75, 12.0, 12.5, 13.0, 13.5,
    14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
    17.0, 17.5, 18.0, 18.5, 19.0, 19.5,
    20.0, 20.5, 21.0, 21.5, 22.0, 22.5,
    23.0, 23.5, 24.0, 24.5, 25.0, 25.5,
    26.0, 26.5, 27.0, 27.5, 28.0, 28.5,
    29.0, 29.5, 30.0, 31.0, 32.0, 33.0,
    34.0, 35.0, 36.0, 37.0, 38.0, 40.0,
    42.0, 44.0, 46.0, 48.0, 50.0, 55.0,
    60.0, 67.5, 75.0, 85.0, 95.0, 107.5,
    120.0, 160.0, 200.0
])

# ============================================================
# 4. Final analysis binning
#
# 34 pT bins x 6 |y| bins
# ============================================================

analysisYBins = array("d", [
    0.0, 0.3, 0.6, 0.9, 1.2, 1.8, 2.4
])

analysisPtBins = array("d", [
    10.0, 10.5, 11.0, 11.5, 12.0, 13.0,
    14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
    20.0, 21.0, 22.0, 23.0, 24.0, 25.0,
    26.0, 27.0, 28.0, 29.0, 30.0, 32.0,
    34.0, 36.0, 38.0, 42.0, 46.0, 50.0,
    60.0, 75.0, 95.0, 120.0, 200.0
])

nFinePtBins = len(finePtBins) - 1
nFineYBins = len(fineYBins) - 1
nAnalysisPtBins = len(analysisPtBins) - 1
nAnalysisYBins = len(analysisYBins) - 1

print("\nBinning check")
print("Fine weight map =", nFinePtBins, "x", nFineYBins)
print("Final acceptance =", nAnalysisPtBins, "x", nAnalysisYBins)

if hRef.GetNbinsX() != nFinePtBins:
    raise RuntimeError(
        "hPromptPtY X bin count does not match finePtBins. "
        "Check the reference histogram binning."
    )

if hRef.GetNbinsY() != nFineYBins:
    raise RuntimeError(
        "hPromptPtY Y bin count does not match fineYBins. "
        "Check the reference histogram binning."
    )


def check_axis_edges(axis, edges, axisName):
    """Check every edge to prevent using a mismatched weight map."""

    for edgeIndex, edgeValue in enumerate(edges):
        if edgeIndex == 0:
            actualValue = axis.GetBinLowEdge(1)
        else:
            actualValue = axis.GetBinUpEdge(edgeIndex)

        if abs(actualValue - edgeValue) > 1.0e-6:
            raise RuntimeError(
                "{} binning mismatch at edge {}: "
                "reference={}, configured={}".format(
                    axisName, edgeIndex, actualValue, edgeValue
                )
            )


check_axis_edges(hRef.GetXaxis(), finePtBins, "pT")
check_axis_edges(hRef.GetYaxis(), fineYBins, "|y|")

# ============================================================
# 5. GEN particle configuration
# ============================================================

genPartsMod = "genParticles"
genPartsType = Handle("std::vector<reco::GenParticle>")

# ============================================================
# 6. First pass: build fine-binned GEN distribution
# ============================================================

hGenPtY = ROOT.TH2D(
    "hGenPtY",
    "GEN J/#psi;p_{T}(J/#psi) [GeV];|y(J/#psi)|",
    nFinePtBins,
    finePtBins,
    nFineYBins,
    fineYBins
)

hGenPtY.Sumw2()

nTotalEvents = 0
nGenJpsi = 0
nGenJpsiInWeightRange = 0

print("\n" + "=" * 70)
print("First pass: building fine-binned hGenPtY")
print("=" * 70)

for fileIndex, inputFile in enumerate(files):

    print(
        "First pass: file {}/{}: {}".format(
            fileIndex + 1, len(files), inputFile
        )
    )

    events = Events(inputFile)

    for event in events:

        nTotalEvents += 1

        event.getByLabel(genPartsMod, genPartsType)
        genParts = genPartsType.product()

        for particle in genParts:

            if abs(particle.pdgId()) != 443:
                continue

            if not particle.isLastCopy():
                continue

            nGenJpsi += 1

            pt = particle.pt()
            rapidity = abs(particle.rapidity())

            if not (
                finePtBins[0] <= pt < finePtBins[-1]
            ):
                continue

            if not (
                fineYBins[0] <= rapidity < fineYBins[-1]
            ):
                continue

            hGenPtY.Fill(pt, rapidity)
            nGenJpsiInWeightRange += 1

# ============================================================
# 7. Construct the fine-binned weight map
#
# weight = (data bin content / data integral)
#          ----------------------------------
#          (MC bin content / MC integral)
# ============================================================

refIntegral = hRef.Integral()
genIntegral = hGenPtY.Integral()

if refIntegral <= 0:
    raise RuntimeError("Reference histogram has zero integral.")

if genIntegral <= 0:
    raise RuntimeError("GEN histogram has zero integral.")

hWeight = hRef.Clone("hPtYWeight")
hWeight.SetDirectory(0)
hWeight.Reset()
hWeight.SetTitle(
    "p_{T}-|y| reweighting factor;p_{T} [GeV];|y|"
)

nEmptyGenWeightBins = 0
nEmptyReferenceWeightBins = 0

for ix in range(1, hRef.GetNbinsX() + 1):

    for iy in range(1, hRef.GetNbinsY() + 1):

        refContent = hRef.GetBinContent(ix, iy)
        genContent = hGenPtY.GetBinContent(ix, iy)

        if refContent <= 0:
            nEmptyReferenceWeightBins += 1

        if genContent <= 0:
            nEmptyGenWeightBins += 1
            weight = 0.0

        else:
            refFraction = refContent / refIntegral
            genFraction = genContent / genIntegral
            weight = refFraction / genFraction

        hWeight.SetBinContent(ix, iy, weight)

print("\nWeight map")
print("Reference integral       =", refIntegral)
print("GEN integral             =", genIntegral)
print("Weight map bins          =",
      hWeight.GetNbinsX(), "x", hWeight.GetNbinsY())
print("Empty GEN weight bins    =", nEmptyGenWeightBins)
print("Empty reference bins     =", nEmptyReferenceWeightBins)
print("Minimum weight           =", hWeight.GetMinimum())
print("Maximum weight           =", hWeight.GetMaximum())

# ============================================================
# 8. Weight lookup
# ============================================================

def getPtYWeight(pt, rapidity):

    absoluteRapidity = abs(rapidity)

    if not (
        finePtBins[0] <= pt < finePtBins[-1]
    ):
        return 0.0

    if not (
        fineYBins[0] <= absoluteRapidity < fineYBins[-1]
    ):
        return 0.0

    xBin = hWeight.GetXaxis().FindFixBin(pt)
    yBin = hWeight.GetYaxis().FindFixBin(absoluteRapidity)

    if xBin < 1 or xBin > hWeight.GetNbinsX():
        return 0.0

    if yBin < 1 or yBin > hWeight.GetNbinsY():
        return 0.0

    weight = hWeight.GetBinContent(xBin, yBin)

    if weight <= 0:
        return 0.0

    return weight

# ============================================================
# 9. QA histograms: pT and |y|
#
# Original QA includes every J/psi in analysis phase space.
# Weighted QA includes only events with a positive weight.
# ============================================================

hPtOriginal = ROOT.TH1D(
    "hPtOriginal",
    "GEN J/#psi p_{T};p_{T} [GeV];Events",
    40, 10.0, 200.0
)

hPtWeighted = ROOT.TH1D(
    "hPtWeighted",
    "Weighted GEN J/#psi p_{T};p_{T} [GeV];Weighted events",
    40, 10.0, 200.0
)

hYOriginal = ROOT.TH1D(
    "hYOriginal",
    "GEN J/#psi |y|;|y|;Events",
    30, 0.0, 2.4
)

hYWeighted = ROOT.TH1D(
    "hYWeighted",
    "Weighted GEN J/#psi |y|;|y|;Weighted events",
    30, 0.0, 2.4
)

# ============================================================
# 10. Final acceptance histograms
#
# All four histograms use analysis binning: 34 x 6.
# ============================================================

def make_analysis_histogram(name, title):
    histogram = ROOT.TH2D(
        name,
        title,
        nAnalysisPtBins,
        analysisPtBins,
        nAnalysisYBins,
        analysisYBins
    )
    histogram.Sumw2()
    return histogram


hDenomOriginal = make_analysis_histogram(
    "hDenomOriginal",
    "Unweighted GEN J/#psi denominator;"
    "p_{T}(J/#psi) [GeV];|y(J/#psi)|"
)

hNumerOriginal = make_analysis_histogram(
    "hNumerOriginal",
    "Unweighted accepted J/#psi numerator;"
    "p_{T}(J/#psi) [GeV];|y(J/#psi)|"
)

hDenomWeighted = make_analysis_histogram(
    "hDenomWeighted",
    "Weighted GEN J/#psi denominator;"
    "p_{T}(J/#psi) [GeV];|y(J/#psi)|"
)

hNumerWeighted = make_analysis_histogram(
    "hNumerWeighted",
    "Weighted accepted J/#psi numerator;"
    "p_{T}(J/#psi) [GeV];|y(J/#psi)|"
)

# ============================================================
# 11. Muon acceptance
# ============================================================

def muon_pass_criteria(muon):

    if abs(muon.eta()) > 2.4:
        return False

    if muon.pt() < 5.0:
        return False

    return True


# ============================================================
# 12. Find final-state muons recursively
#
# Preserve the recursive approach from the original code.
# This allows muons to be found below intermediate daughters.
# ============================================================

def find_final_muons(particle):

    muons = []

    if abs(particle.pdgId()) == 13:

        if particle.numberOfDaughters() == 0:
            muons.append(particle)

        else:
            for daughterIndex in range(
                particle.numberOfDaughters()
            ):
                daughter = particle.daughter(daughterIndex)

                if daughter is not None:
                    muons.extend(find_final_muons(daughter))

        return muons

    for daughterIndex in range(
        particle.numberOfDaughters()
    ):
        daughter = particle.daughter(daughterIndex)

        if daughter is not None:
            muons.extend(find_final_muons(daughter))

    return muons


# ============================================================
# 13. Second pass: calculate acceptance
#
# IMPORTANT:
# - Original denominator is filled independently of weight.
# - Weighted denominator is filled only when weight > 0.
# - Original numerator is also independent of weight.
# - Weighted numerator requires a positive weight.
# ============================================================

nAnalysisJpsi = 0
nWeightedJpsi = 0
nAcceptedOriginal = 0
nAcceptedWeighted = 0

sumWeights = 0.0
sumWeightsAccepted = 0.0

print("\n" + "=" * 70)
print("Second pass: calculating acceptance")
print("=" * 70)

for fileIndex, inputFile in enumerate(files):

    print(
        "Second pass: file {}/{}: {}".format(
            fileIndex + 1, len(files), inputFile
        )
    )

    events = Events(inputFile)

    for event in events:

        event.getByLabel(genPartsMod, genPartsType)
        genParts = genPartsType.product()

        for particle in genParts:

            if abs(particle.pdgId()) != 443:
                continue

            if not particle.isLastCopy():
                continue

            pt = particle.pt()
            rapidity = abs(particle.rapidity())

            # ------------------------------------------------
            # Final analysis phase space
            # ------------------------------------------------

            if not (
                analysisPtBins[0] <= pt < analysisPtBins[-1]
            ):
                continue

            if not (
                analysisYBins[0] <= rapidity < analysisYBins[-1]
            ):
                continue

            nAnalysisJpsi += 1

            # ------------------------------------------------
            # Original, unweighted denominator
            # ------------------------------------------------

            hDenomOriginal.Fill(pt, rapidity)
            hPtOriginal.Fill(pt)
            hYOriginal.Fill(rapidity)

            # ------------------------------------------------
            # Determine whether the decay passes acceptance.
            # Do this independently of the reweighting value.
            # ------------------------------------------------

            muons = find_final_muons(particle)

            accepted = False

            if len(muons) == 2:
                accepted = (
                    muon_pass_criteria(muons[0])
                    and muon_pass_criteria(muons[1])
                )

            # ------------------------------------------------
            # Original, unweighted numerator
            # ------------------------------------------------

            if accepted:
                hNumerOriginal.Fill(pt, rapidity)
                nAcceptedOriginal += 1

            # ------------------------------------------------
            # Weighted denominator and QA
            # ------------------------------------------------

            weight = getPtYWeight(pt, rapidity)

            if weight <= 0:
                continue

            nWeightedJpsi += 1
            sumWeights += weight

            hDenomWeighted.Fill(pt, rapidity, weight)
            hPtWeighted.Fill(pt, weight)
            hYWeighted.Fill(rapidity, weight)

            # ------------------------------------------------
            # Weighted numerator
            # ------------------------------------------------

            if accepted:
                hNumerWeighted.Fill(pt, rapidity, weight)
                nAcceptedWeighted += 1
                sumWeightsAccepted += weight

# ============================================================
# 14. Calculate bin-by-bin acceptance
# ============================================================

hAcceptanceOriginal = hNumerOriginal.Clone(
    "hAcceptanceOriginal"
)
hAcceptanceOriginal.SetDirectory(0)
hAcceptanceOriginal.SetTitle(
    "Unweighted J/#psi acceptance;"
    "p_{T}(J/#psi) [GeV];|y(J/#psi)|"
)
hAcceptanceOriginal.Reset()
hAcceptanceOriginal.Divide(
    hNumerOriginal,
    hDenomOriginal,
    1.0,
    1.0,
    ""
)

hAcceptanceWeighted = hNumerWeighted.Clone(
    "hAcceptanceWeighted"
)
hAcceptanceWeighted.SetDirectory(0)
hAcceptanceWeighted.SetTitle(
    "Weighted J/#psi acceptance;"
    "p_{T}(J/#psi) [GeV];|y(J/#psi)|"
)
hAcceptanceWeighted.Reset()
hAcceptanceWeighted.Divide(
    hNumerWeighted,
    hDenomWeighted,
    1.0,
    1.0,
    ""
)

# ============================================================
# 15. Integrated acceptance and uncertainties
#
# Sumw2 errors are used for weighted histograms.
# The ratio uncertainty below uses independent numerator and
# denominator error propagation, matching the original method.
# ============================================================

def get_integrated_acceptance(denominatorHistogram,
                              numeratorHistogram):

    denominator = denominatorHistogram.Integral()
    numerator = numeratorHistogram.Integral()

    denominatorErrorSquared = 0.0
    numeratorErrorSquared = 0.0

    for xBin in range(
        1, denominatorHistogram.GetNbinsX() + 1
    ):
        for yBin in range(
            1, denominatorHistogram.GetNbinsY() + 1
        ):
            denominatorError = denominatorHistogram.GetBinError(
                xBin, yBin
            )
            numeratorError = numeratorHistogram.GetBinError(
                xBin, yBin
            )

            denominatorErrorSquared += denominatorError ** 2
            numeratorErrorSquared += numeratorError ** 2

    denominatorError = math.sqrt(denominatorErrorSquared)
    numeratorError = math.sqrt(numeratorErrorSquared)

    if denominator <= 0:
        return (
            denominator, numerator,
            denominatorError, numeratorError,
            0.0, 0.0
        )

    acceptance = numerator / denominator

    if numerator > 0:
        acceptanceError = acceptance * math.sqrt(
            (numeratorError / numerator) ** 2
            +
            (denominatorError / denominator) ** 2
        )
    else:
        acceptanceError = 0.0

    return (
        denominator,
        numerator,
        denominatorError,
        numeratorError,
        acceptance,
        acceptanceError
    )


(
    totalDenomOriginal,
    totalNumerOriginal,
    errorDenomOriginal,
    errorNumerOriginal,
    integratedAcceptanceOriginal,
    integratedAcceptanceOriginalError
) = get_integrated_acceptance(
    hDenomOriginal,
    hNumerOriginal
)

(
    totalDenomWeighted,
    totalNumerWeighted,
    errorDenomWeighted,
    errorNumerWeighted,
    integratedAcceptanceWeighted,
    integratedAcceptanceWeightedError
) = get_integrated_acceptance(
    hDenomWeighted,
    hNumerWeighted
)

# ============================================================
# 16. Consistency checks
# ============================================================

print("\n" + "=" * 70)
print("Counter / histogram consistency checks")
print("=" * 70)

print("Analysis J/psi counter       =", nAnalysisJpsi)
print("Original denominator integral =", totalDenomOriginal)

print("Original accepted counter     =", nAcceptedOriginal)
print("Original numerator integral   =", totalNumerOriginal)

print("Positive-weight counter       =", nWeightedJpsi)
print("Weighted denominator integral =", totalDenomWeighted)
print("Sum of weights                =", sumWeights)

print("Weighted accepted counter     =", nAcceptedWeighted)
print("Weighted numerator integral   =", totalNumerWeighted)
print("Accepted sum of weights       =", sumWeightsAccepted)

print("Check original denominator:",
      abs(nAnalysisJpsi - totalDenomOriginal))
print("Check original numerator:",
      abs(nAcceptedOriginal - totalNumerOriginal))
print("Check weighted denominator:",
      abs(sumWeights - totalDenomWeighted))
print("Check weighted numerator:",
      abs(sumWeightsAccepted - totalNumerWeighted))

print("=" * 70)

# ============================================================
# 17. Print final summary
# ============================================================

print("\n" + "=" * 70)
print("Final acceptance summary")
print("=" * 70)

print("Total input events              =", nTotalEvents)
print("All last-copy GEN J/psi         =", nGenJpsi)
print("GEN J/psi in fine weight range  =", nGenJpsiInWeightRange)
print("GEN J/psi in analysis range     =", nAnalysisJpsi)

print()
print("Original denominator = {:.8f} +/- {:.8f}".format(
    totalDenomOriginal, errorDenomOriginal
))
print("Original numerator   = {:.8f} +/- {:.8f}".format(
    totalNumerOriginal, errorNumerOriginal
))
print("Original acceptance  = {:.8f} +/- {:.8f}".format(
    integratedAcceptanceOriginal,
    integratedAcceptanceOriginalError
))

print()
print("Weighted denominator = {:.8f} +/- {:.8f}".format(
    totalDenomWeighted, errorDenomWeighted
))
print("Weighted numerator   = {:.8f} +/- {:.8f}".format(
    totalNumerWeighted, errorNumerWeighted
))
print("Weighted acceptance  = {:.8f} +/- {:.8f}".format(
    integratedAcceptanceWeighted,
    integratedAcceptanceWeightedError
))

print()
print("Sum weights          = {:.8f}".format(sumWeights))
print("Accepted sum weights = {:.8f}".format(sumWeightsAccepted))

print("Acceptance difference = {:.8f}".format(
    integratedAcceptanceWeighted - integratedAcceptanceOriginal
))
print("=" * 70)

# ============================================================
# 18. Write detailed TXT results
# ============================================================

with open(outputTextFileName, "w") as outputTextFile:

    outputTextFile.write("=" * 60 + "\n")
    outputTextFile.write(
        "Numerator / Denominator  (unweighted and weighted)\n"
    )
    outputTextFile.write("=" * 60 + "\n")

    outputTextFile.write(
        "Number of input files       = {}\n".format(len(files))
    )
    outputTextFile.write(
        "Total events                 = {}\n".format(nTotalEvents)
    )
    outputTextFile.write(
        "Total GEN J/psi              = {}\n".format(nGenJpsi)
    )
    outputTextFile.write(
        "Selected J/psi               = {}\n".format(nAnalysisJpsi)
    )
    outputTextFile.write(
        "Accepted J/psi               = {}\n".format(nAcceptedOriginal)
    )
    outputTextFile.write(
        "Weighted denominator         = {:.10f}\n".format(
            hDenomWeighted.Integral()
        )
    )
    outputTextFile.write(
        "Weighted numerator           = {:.10f}\n".format(
            hNumerWeighted.Integral()
        )
    )
    outputTextFile.write(
        "Sum of weights               = {:.10f}\n".format(
            sumWeights
        )
    )
    outputTextFile.write(
        "Integrated acceptance (unw)  = {:.10f}\n".format(
            integratedAcceptanceOriginal
        )
    )
    outputTextFile.write(
        "Integrated acceptance (w)   = {:.10f}\n".format(
            integratedAcceptanceWeighted
        )
    )

    outputTextFile.write("\n")
    outputTextFile.write("=" * 60 + "\n")
    outputTextFile.write(
        "Bin-by-bin numerator and denominator\n"
    )
    outputTextFile.write("=" * 60 + "\n")

    outputTextFile.write(
        "{:<12s}{:<13s}{:<12s}{:<12s}"
        "{:<20s}{:<20s}{:<20s}"
        "{:<20s}{:<20s}{:<20s}\n".format(
            "pT_low", "pT_high", "y_low", "y_high",
            "Denominator", "Numerator", "Acceptance",
            "Den_w", "Num_w", "Acc_w"
        )
    )
    outputTextFile.write("-" * 160 + "\n")

    for ptBin in range(1, len(analysisPtBins)):
        for yBin in range(1, len(analysisYBins)):

            ptLow = analysisPtBins[ptBin - 1]
            ptHigh = analysisPtBins[ptBin]
            yLow = analysisYBins[yBin - 1]
            yHigh = analysisYBins[yBin]

            denominator = hDenomOriginal.GetBinContent(
                ptBin, yBin
            )
            numerator = hNumerOriginal.GetBinContent(
                ptBin, yBin
            )
            denominatorWeighted = hDenomWeighted.GetBinContent(
                ptBin, yBin
            )
            numeratorWeighted = hNumerWeighted.GetBinContent(
                ptBin, yBin
            )

            acceptance = (
                numerator / denominator
                if denominator > 0 else 0.0
            )
            acceptanceWeighted = (
                numeratorWeighted / denominatorWeighted
                if denominatorWeighted > 0 else 0.0
            )

            outputTextFile.write(
                "{:<12.3f}{:<13.3f}{:<12.3f}{:<12.3f}"
                "{:<20.8f}{:<20.8f}{:<20.8f}"
                "{:<20.8f}{:<20.8f}{:<20.8f}\n".format(
                    ptLow, ptHigh, yLow, yHigh,
                    denominator, numerator, acceptance,
                    denominatorWeighted, numeratorWeighted,
                    acceptanceWeighted
                )
            )

# ============================================================
# 19. Draw pT and |y| QA comparison plots
# ============================================================

def draw_compare(histogramOriginal, histogramWeighted, canvasName):

    canvas = ROOT.TCanvas(canvasName, canvasName, 900, 700)

    histogramOriginal.SetLineColor(ROOT.kBlack)
    histogramOriginal.SetLineWidth(2)

    histogramWeighted.SetLineColor(ROOT.kRed)
    histogramWeighted.SetLineWidth(2)

    maximumValue = max(
        histogramOriginal.GetMaximum(),
        histogramWeighted.GetMaximum()
    )

    if maximumValue > 0:
        histogramOriginal.SetMaximum(maximumValue * 1.35)

    histogramOriginal.Draw("HIST")
    histogramWeighted.Draw("HIST SAME")

    legend = ROOT.TLegend(0.18, 0.75, 0.48, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(histogramOriginal, "Particle gun", "l")
    legend.AddEntry(
        histogramWeighted,
        "After p_{T}-|y| weighting",
        "l"
    )
    legend.Draw()

    canvas.Update()
    canvas.SaveAs(canvasName + ".png")

    return canvas


cPtCompare = draw_compare(
    hPtOriginal,
    hPtWeighted,
    "JpsiPt_before_after_v3"
)

cYCompare = draw_compare(
    hYOriginal,
    hYWeighted,
    "JpsiY_before_after_v3"
)

# ============================================================
# 20. Draw fine-binned weight map
# ============================================================

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat("4.3f")

cWeight = ROOT.TCanvas(
    "cWeight",
    "pT-y weight map",
    1000,
    800
)

hWeight.SetTitle(
    "p_{T}-|y| reweighting factor;p_{T} [GeV];|y|"
)
hWeight.SetMarkerSize(0.7)
hWeight.SetMinimum(0.0)
hWeight.Draw("COLZ TEXT")

cWeight.Update()
cWeight.SaveAs("Jpsi_pT_y_weight_v3.png")

# ============================================================
# 21. Transpose acceptance histograms
#
# Original: X = pT, Y = |y|
# Transposed: X = |y|, Y = pT
# ============================================================

hAcceptanceWeighted_Transpose = ROOT.TH2D(
    "hAcceptanceWeighted_Transpose",
    "Weighted J/#psi acceptance;|y|;p_{T} [GeV]",
    nAnalysisYBins,
    analysisYBins,
    nAnalysisPtBins,
    analysisPtBins
)

hAcceptanceOriginal_Transpose = ROOT.TH2D(
    "hAcceptanceOriginal_Transpose",
    "Unweighted J/#psi acceptance;|y|;p_{T} [GeV]",
    nAnalysisYBins,
    analysisYBins,
    nAnalysisPtBins,
    analysisPtBins
)

hAcceptanceWeighted_Transpose.Sumw2()
hAcceptanceOriginal_Transpose.Sumw2()

for xBin in range(1, nAnalysisPtBins + 1):

    for yBin in range(1, nAnalysisYBins + 1):

        weightedValue = hAcceptanceWeighted.GetBinContent(
            xBin, yBin
        )
        weightedError = hAcceptanceWeighted.GetBinError(
            xBin, yBin
        )

        originalValue = hAcceptanceOriginal.GetBinContent(
            xBin, yBin
        )
        originalError = hAcceptanceOriginal.GetBinError(
            xBin, yBin
        )

        hAcceptanceWeighted_Transpose.SetBinContent(
            yBin, xBin, weightedValue
        )
        hAcceptanceWeighted_Transpose.SetBinError(
            yBin, xBin, weightedError
        )

        hAcceptanceOriginal_Transpose.SetBinContent(
            yBin, xBin, originalValue
        )
        hAcceptanceOriginal_Transpose.SetBinError(
            yBin, xBin, originalError
        )

# ============================================================
# 22. Draw weighted acceptance map
# ============================================================

ROOT.gStyle.SetPaintTextFormat("4.2f")

cAcceptanceWeighted = ROOT.TCanvas(
    "cAcceptanceWeighted",
    "Weighted acceptance map",
    1000,
    800
)

cAcceptanceWeighted.SetLogy()

hAcceptanceWeighted_Transpose.SetTitle(
    "Weighted J/#psi acceptance;|y|;p_{T} [GeV]"
)
hAcceptanceWeighted_Transpose.SetMarkerSize(0.7)
hAcceptanceWeighted_Transpose.SetMinimum(0.0)
hAcceptanceWeighted_Transpose.SetMaximum(1.0)
hAcceptanceWeighted_Transpose.GetXaxis().SetTitle("|y|")
hAcceptanceWeighted_Transpose.GetYaxis().SetTitle("p_{T} [GeV]")
hAcceptanceWeighted_Transpose.GetXaxis().SetRangeUser(0.0, 2.4)
hAcceptanceWeighted_Transpose.GetYaxis().SetRangeUser(10.0, 200.0)

hAcceptanceWeighted_Transpose.Draw("COLZ")
cAcceptanceWeighted.Update()
cAcceptanceWeighted.SaveAs("Jpsi_acceptance_weighted_v3.png")

# ============================================================
# 23. Draw original acceptance map
# ============================================================

cAcceptanceOriginal = ROOT.TCanvas(
    "cAcceptanceOriginal",
    "Original acceptance map",
    1000,
    800
)

cAcceptanceOriginal.SetLogy()

hAcceptanceOriginal_Transpose.SetTitle(
    "Unweighted J/#psi acceptance;|y|;p_{T} [GeV]"
)
hAcceptanceOriginal_Transpose.SetMarkerSize(0.7)
hAcceptanceOriginal_Transpose.SetMinimum(0.0)
hAcceptanceOriginal_Transpose.SetMaximum(1.0)
hAcceptanceOriginal_Transpose.GetXaxis().SetTitle("|y|")
hAcceptanceOriginal_Transpose.GetYaxis().SetTitle("p_{T} [GeV]")
hAcceptanceOriginal_Transpose.GetXaxis().SetRangeUser(0.0, 2.4)
hAcceptanceOriginal_Transpose.GetYaxis().SetRangeUser(10.0, 200.0)

hAcceptanceOriginal_Transpose.Draw("COLZ")
cAcceptanceOriginal.Update()
cAcceptanceOriginal.SaveAs("Jpsi_acceptance_original_v3.png")

# ============================================================
# 24. Save all important histograms to ROOT
# ============================================================

outputFile = ROOT.TFile(outputRootFileName, "RECREATE")

if not outputFile or outputFile.IsZombie():
    raise RuntimeError("Cannot create output ROOT file: " + outputRootFileName)

histogramsToWrite = [
    hRef,
    hGenPtY,
    hWeight,

    hPtOriginal,
    hPtWeighted,
    hYOriginal,
    hYWeighted,

    hDenomOriginal,
    hNumerOriginal,
    hDenomWeighted,
    hNumerWeighted,

    hAcceptanceOriginal,
    hAcceptanceWeighted,

    hAcceptanceOriginal_Transpose,
    hAcceptanceWeighted_Transpose
]

for histogram in histogramsToWrite:
    histogram.Write()

outputFile.Close()
fRef.Close()

# ============================================================
# 25. Finished
# ============================================================

print("\n" + "=" * 70)
print("Finished")
print("=" * 70)
print("Events                         =", nTotalEvents)
print("All last-copy GEN J/psi        =", nGenJpsi)
print("GEN J/psi in analysis range    =", nAnalysisJpsi)
print("Original accepted J/psi        =", nAcceptedOriginal)
print("Positive-weight J/psi          =", nWeightedJpsi)
print("Weighted accepted J/psi        =", nAcceptedWeighted)
print("Original acceptance            =", integratedAcceptanceOriginal)
print("Weighted acceptance            =", integratedAcceptanceWeighted)
print()
print("ROOT output:")
print(" ", outputRootFileName)
print()
print("TXT output:")
print(" ", outputTextFileName)
print()
print("Plots:")
print("  JpsiPt_before_after_v3.png")
print("  JpsiY_before_after_v3.png")
print("  Jpsi_pT_y_weight_v3.png")
print("  Jpsi_acceptance_weighted_v3.png")
print("  Jpsi_acceptance_original_v3.png")
print("=" * 70)

