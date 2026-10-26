#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ROOT
import glob
from array import array
from DataFormats.FWLite import Events, Handle

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# ============================================================
# Input files
# ============================================================

refFile = "fitOutput_final.root"

files = glob.glob(
    "/eos/user/z/zhipeng/MC_Particle_Gun_For_Acc/results/step0_*.root"
)

print("============================================================")
print("Input files")
print("============================================================")
print("Number of files =", len(files))
print("Reference file  =", refFile)
print("============================================================")


# ============================================================
# Read reference histogram
# This histogram defines the FINE binning for the weight map
# ============================================================

fRef = ROOT.TFile.Open(refFile)

if not fRef or fRef.IsZombie():
    raise RuntimeError("Cannot open reference file: " + refFile)

hRef = fRef.Get("hPromptPtY")

if not hRef:
    raise RuntimeError(
        "Cannot find hPromptPtY in " + refFile
    )

print()
print("============================================================")
print("Reference histogram")
print("============================================================")
print("Name      =", hRef.GetName())
print("Integral  =", hRef.Integral())
print("Nbins X   =", hRef.GetNbinsX())
print("Nbins Y   =", hRef.GetNbinsY())
print("X range   =",
      hRef.GetXaxis().GetXmin(),
      "to",
      hRef.GetXaxis().GetXmax())
print("Y range   =",
      hRef.GetYaxis().GetXmin(),
      "to",
      hRef.GetYaxis().GetXmax())
print("============================================================")


# ============================================================
# Fine binning
#
# These bins come from hPromptPtY.
# They are used ONLY for the pT-y weight map.
#
# 68 pT bins × 12 |y| bins
# ============================================================

fineYBins = array("d", [
    0.0,
    0.15,
    0.3,
    0.45,
    0.6,
    0.75,
    0.9,
    1.05,
    1.2,
    1.5,
    1.8,
    2.1,
    2.4
])

finePtBins = array("d", [
    10.0,
    10.25,
    10.5,
    10.75,
    11.0,
    11.25,
    11.5,
    11.75,
    12.0,
    12.5,
    13.0,
    13.5,
    14.0,
    14.5,
    15.0,
    15.5,
    16.0,
    16.5,
    17.0,
    17.5,
    18.0,
    18.5,
    19.0,
    19.5,
    20.0,
    20.5,
    21.0,
    21.5,
    22.0,
    22.5,
    23.0,
    23.5,
    24.0,
    24.5,
    25.0,
    25.5,
    26.0,
    26.5,
    27.0,
    27.5,
    28.0,
    28.5,
    29.0,
    29.5,
    30.0,
    31.0,
    32.0,
    33.0,
    34.0,
    35.0,
    36.0,
    37.0,
    38.0,
    40.0,
    42.0,
    44.0,
    46.0,
    48.0,
    50.0,
    55.0,
    60.0,
    67.5,
    75.0,
    85.0,
    95.0,
    107.5,
    120.0,
    160.0,
    200.0
])


# ============================================================
# Final analysis binning
#
# These bins are used ONLY for the final acceptance.
#
# 34 pT bins × 6 |y| bins
# ============================================================

analysisYBins = array("d", [
    0.0,
    0.3,
    0.6,
    0.9,
    1.2,
    1.8,
    2.4
])

analysisPtBins = array("d", [
    10.0,
    10.5,
    11.0,
    11.5,
    12.0,
    13.0,
    14.0,
    15.0,
    16.0,
    17.0,
    18.0,
    19.0,
    20.0,
    21.0,
    22.0,
    23.0,
    24.0,
    25.0,
    26.0,
    27.0,
    28.0,
    29.0,
    30.0,
    32.0,
    34.0,
    36.0,
    38.0,
    42.0,
    46.0,
    50.0,
    60.0,
    75.0,
    95.0,
    120.0,
    200.0
])


# ============================================================
# Check binning
# ============================================================

print()
print("============================================================")
print("Binning check")
print("============================================================")
print("Fine pT bins     =", len(finePtBins) - 1)
print("Fine |y| bins    =", len(fineYBins) - 1)
print("Fine 2D bins     =",
      (len(finePtBins) - 1) * (len(fineYBins) - 1))

print()

print("Analysis pT bins =", len(analysisPtBins) - 1)
print("Analysis |y| bins=", len(analysisYBins) - 1)
print("Analysis 2D bins =",
      (len(analysisPtBins) - 1) * (len(analysisYBins) - 1))

print("============================================================")


# ============================================================
# GEN particle information
# ============================================================

genPartsMod = "genParticles"
genPartsType = Handle("vector<reco::GenParticle>")


# ============================================================
# First pass
#
# Build MC GEN pT-y distribution using the FINE bins.
#
# This histogram must have exactly the same binning as hRef.
# ============================================================

hGenPtY = ROOT.TH2D(
    "hGenPtY",
    "GEN J/#psi;p_{T}(J/#psi) [GeV];|y(J/#psi)|",
    len(finePtBins) - 1,
    finePtBins,
    len(fineYBins) - 1,
    fineYBins
)

nTotalEvents = 0
nGenJpsi = 0


print()
print("============================================================")
print("First pass: building hGenPtY")
print("============================================================")


for fileIndex, inputFile in enumerate(files):

    print(
        "Processing file",
        fileIndex + 1,
        "/",
        len(files),
        ":",
        inputFile
    )

    events = Events(inputFile)

    for event in events:

        nTotalEvents += 1

        event.getByLabel(
            genPartsMod,
            genPartsType
        )

        genParts = genPartsType.product()

        for p in genParts:

            if abs(p.pdgId()) != 443:
                continue

            if not p.isLastCopy():
                continue

            pt = p.pt()
            y = abs(p.rapidity())

            if not (
                hRef.GetXaxis().GetXmin()
                <= pt
                <
                hRef.GetXaxis().GetXmax()
            ):
                continue

            if not (
                hRef.GetYaxis().GetXmin()
                <= y
                <
                hRef.GetYaxis().GetXmax()
            ):
                continue

            hGenPtY.Fill(pt, y)

            nGenJpsi += 1


# ============================================================
# Build weight map
#
# weight = Data fraction / MC fraction
#
# The weight map remains FINE BINNED: 68 × 12
# ============================================================

refIntegral = hRef.Integral()
genIntegral = hGenPtY.Integral()

if refIntegral <= 0:
    raise RuntimeError(
        "Reference histogram integral is zero."
    )

if genIntegral <= 0:
    raise RuntimeError(
        "GEN histogram integral is zero."
    )


hWeight = hRef.Clone("hPtYWeight")
hWeight.Reset()


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

        hWeight.SetBinContent(
            ix,
            iy,
            weight
        )


# ============================================================
# Weight function
# ============================================================

def getPtYWeight(pt, y):

    y = abs(y)

    if not (
        hWeight.GetXaxis().GetXmin()
        <= pt
        <
        hWeight.GetXaxis().GetXmax()
    ):
        return 0.0

    if not (
        hWeight.GetYaxis().GetXmin()
        <= y
        <
        hWeight.GetYaxis().GetXmax()
    ):
        return 0.0

    ix = hWeight.GetXaxis().FindFixBin(pt)
    iy = hWeight.GetYaxis().FindFixBin(y)

    weight = hWeight.GetBinContent(ix, iy)

    if weight <= 0:
        return 0.0

    return weight


# ============================================================
# Print weight information
# ============================================================

print()
print("============================================================")
print("Weight map")
print("============================================================")
print("Reference integral =", refIntegral)
print("GEN integral       =", genIntegral)
print("Weight map bins X  =", hWeight.GetNbinsX())
print("Weight map bins Y  =", hWeight.GetNbinsY())
print("============================================================")


# ============================================================
# QA histograms
#
# These are not used for the final acceptance.
# ============================================================

hPtOriginal = ROOT.TH1D(
    "hPtOriginal",
    "GEN J/#psi p_{T};p_{T} [GeV];Events",
    40,
    10.0,
    200.0
)

hPtWeighted = ROOT.TH1D(
    "hPtWeighted",
    "Weighted GEN J/#psi p_{T};p_{T} [GeV];Weighted Events",
    40,
    10.0,
    200.0
)

hYOriginal = ROOT.TH1D(
    "hYOriginal",
    "GEN J/#psi |y|;|y|;Events",
    30,
    0.0,
    2.4
)

hYWeighted = ROOT.TH1D(
    "hYWeighted",
    "Weighted GEN J/#psi |y|;|y|;Weighted Events",
    30,
    0.0,
    2.4
)


# ============================================================
# Final acceptance histograms
#
# IMPORTANT:
#
# These are ANALYSIS BINNED: 34 × 6
#
# The weight itself is still taken from the 68 × 12
# fine-bin weight map.
# ============================================================

hDenomOriginal = ROOT.TH2D(
    "hDenomOriginal",
    "Unweighted GEN J/#psi;p_{T}(J/#psi) [GeV];|y(J/#psi)|",
    len(analysisPtBins) - 1,
    analysisPtBins,
    len(analysisYBins) - 1,
    analysisYBins
)

hNumerOriginal = ROOT.TH2D(
    "hNumerOriginal",
    "Unweighted accepted J/#psi;p_{T}(J/#psi) [GeV];|y(J/#psi)|",
    len(analysisPtBins) - 1,
    analysisPtBins,
    len(analysisYBins) - 1,
    analysisYBins
)

hDenomWeighted = ROOT.TH2D(
    "hDenomWeighted",
    "Weighted GEN J/#psi;p_{T}(J/#psi) [GeV];|y(J/#psi)|",
    len(analysisPtBins) - 1,
    analysisPtBins,
    len(analysisYBins) - 1,
    analysisYBins
)

hNumerWeighted = ROOT.TH2D(
    "hNumerWeighted",
    "Weighted accepted J/#psi;p_{T}(J/#psi) [GeV];|y(J/#psi)|",
    len(analysisPtBins) - 1,
    analysisPtBins,
    len(analysisYBins) - 1,
    analysisYBins
)


print()
print("============================================================")
print("Acceptance histogram")
print("============================================================")
print("Denom X =", hDenomWeighted.GetNbinsX())
print("Denom Y =", hDenomWeighted.GetNbinsY())
print("Numer X =", hNumerWeighted.GetNbinsX())
print("Numer Y =", hNumerWeighted.GetNbinsY())
print("============================================================")


# ============================================================
# Second pass
#
# Fill final acceptance using ANALYSIS bins.
#
# Weight is obtained from the FINE 68 × 12 map.
# ============================================================

nSelected = 0
nAccepted = 0

sumWeights = 0.0
sumWeightsAccepted = 0.0


print()
print("============================================================")
print("Second pass: calculating acceptance")
print("============================================================")


for fileIndex, inputFile in enumerate(files):

    print(
        "Processing file",
        fileIndex + 1,
        "/",
        len(files),
        ":",
        inputFile
    )

    events = Events(inputFile)

    for event in events:

        event.getByLabel(
            genPartsMod,
            genPartsType
        )

        genParts = genPartsType.product()

        for p in genParts:

            if abs(p.pdgId()) != 443:
                continue

            if not p.isLastCopy():
                continue

            pt = p.pt()
            y = abs(p.rapidity())

            # ------------------------------------------------
            # Analysis phase space
            # ------------------------------------------------

            if not (
                analysisPtBins[0]
                <= pt
                <
                analysisPtBins[-1]
            ):
                continue

            if not (
                analysisYBins[0]
                <= y
                <
                analysisYBins[-1]
            ):
                continue

            # ------------------------------------------------
            # Get weight from FINE bin
            # ------------------------------------------------

            weight = getPtYWeight(pt, y)

            if weight <= 0:
                continue

            nSelected += 1
            sumWeights += weight

            # ------------------------------------------------
            # QA
            # ------------------------------------------------

            hPtOriginal.Fill(pt)
            hPtWeighted.Fill(pt, weight)

            hYOriginal.Fill(y)
            hYWeighted.Fill(y, weight)

            # ------------------------------------------------
            # Final acceptance denominator
            #
            # NOTE:
            # This histogram uses 34 × 6 analysis bins.
            # ------------------------------------------------

            hDenomOriginal.Fill(
                pt,
                y
            )

            hDenomWeighted.Fill(
                pt,
                y,
                weight
            )

            # ------------------------------------------------
            # Find J/psi decay muons
            # ------------------------------------------------

            muons = []

            for i in range(p.numberOfDaughters()):

                daughter = p.daughter(i)

                if daughter is None:
                    continue

                if abs(daughter.pdgId()) != 13:
                    continue

                muons.append(daughter)

            # ------------------------------------------------
            # Require exactly two muons
            # ------------------------------------------------

            if len(muons) != 2:
                continue

            # ------------------------------------------------
            # Muon acceptance
            # ------------------------------------------------

            accepted = True

            for muon in muons:

                if abs(muon.eta()) > 2.4:
                    accepted = False
                    break

                if muon.pt() < 5.0:
                    accepted = False
                    break

            if not accepted:
                continue

            # ------------------------------------------------
            # Final acceptance numerator
            # ------------------------------------------------

            hNumerOriginal.Fill(
                pt,
                y
            )

            hNumerWeighted.Fill(
                pt,
                y,
                weight
            )

            nAccepted += 1

            sumWeightsAccepted += weight


# ============================================================
# Calculate acceptance
# ============================================================

hAcceptanceOriginal = hNumerOriginal.Clone(
    "hAcceptanceOriginal"
)

hAcceptanceOriginal.SetTitle(
    "Unweighted acceptance;p_{T}(J/#psi) [GeV];|y(J/#psi)|"
)

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

hAcceptanceWeighted.SetTitle(
    "Weighted acceptance;p_{T}(J/#psi) [GeV];|y(J/#psi)|"
)

hAcceptanceWeighted.Divide(
    hNumerWeighted,
    hDenomWeighted,
    1.0,
    1.0,
    ""
)


# ============================================================
# Integrated acceptance
# ============================================================

totalDenomOriginal = hDenomOriginal.Integral()
totalNumerOriginal = hNumerOriginal.Integral()

totalDenomWeighted = hDenomWeighted.Integral()
totalNumerWeighted = hNumerWeighted.Integral()


if totalDenomOriginal > 0:

    integratedAcceptanceOriginal = (
        totalNumerOriginal
        /
        totalDenomOriginal
    )

else:

    integratedAcceptanceOriginal = 0.0


if totalDenomWeighted > 0:

    integratedAcceptanceWeighted = (
        totalNumerWeighted
        /
        totalDenomWeighted
    )

else:

    integratedAcceptanceWeighted = 0.0


# ============================================================
# Print summary
# ============================================================

print()
print("============================================================")
print("Final summary")
print("============================================================")

print("Total events        =", nTotalEvents)
print("Total GEN J/psi     =", nGenJpsi)

print()
print("Selected J/psi      =", nSelected)
print("Accepted J/psi      =", nAccepted)

print()
print("Sum weights         =", sumWeights)
print("Accepted sum weights=",
      sumWeightsAccepted)

print()
print("Original denominator =", totalDenomOriginal)
print("Original numerator   =", totalNumerOriginal)
print("Original acceptance  =",
      integratedAcceptanceOriginal)

print()
print("Weighted denominator =", totalDenomWeighted)
print("Weighted numerator   =", totalNumerWeighted)
print("Weighted acceptance  =",
      integratedAcceptanceWeighted)

print()
print("Acceptance difference =",
      integratedAcceptanceWeighted
      -
      integratedAcceptanceOriginal)

print("============================================================")


# ============================================================
# Create transposed acceptance histogram
#
# Original acceptance:
#   X = pT
#   Y = |y|
#
# Transposed:
#   X = |y|
#   Y = pT
#
# Still uses FINAL ANALYSIS bins.
# ============================================================

hAcceptanceWeighted_Transpose = ROOT.TH2D(
    "hAcceptanceWeighted_Transpose",
    "Weighted acceptance;|y(J/#psi)|;p_{T}(J/#psi) [GeV]",
    len(analysisYBins) - 1,
    analysisYBins,
    len(analysisPtBins) - 1,
    analysisPtBins
)


hAcceptanceOriginal_Transpose = ROOT.TH2D(
    "hAcceptanceOriginal_Transpose",
    "Unweighted acceptance;|y(J/#psi)|;p_{T}(J/#psi) [GeV]",
    len(analysisYBins) - 1,
    analysisYBins,
    len(analysisPtBins) - 1,
    analysisPtBins
)


# ============================================================
# Fill transposed histograms
# ============================================================

for ix in range(
    1,
    hAcceptanceWeighted.GetNbinsX() + 1
):

    for iy in range(
        1,
        hAcceptanceWeighted.GetNbinsY() + 1
    ):

        valueWeighted = (
            hAcceptanceWeighted.GetBinContent(
                ix,
                iy
            )
        )

        valueOriginal = (
            hAcceptanceOriginal.GetBinContent(
                ix,
                iy
            )
        )

        hAcceptanceWeighted_Transpose.SetBinContent(
            iy,
            ix,
            valueWeighted
        )

        hAcceptanceOriginal_Transpose.SetBinContent(
            iy,
            ix,
            valueOriginal
        )


# ============================================================
# Output file
# ============================================================

outputFile = ROOT.TFile(
    "Acceptance_reweighted_v2.root",
    "RECREATE"
)


# Reference / weight information
hRef.Write()
hGenPtY.Write()
hWeight.Write()

# QA
hPtOriginal.Write()
hPtWeighted.Write()
hYOriginal.Write()
hYWeighted.Write()

# Acceptance
hDenomOriginal.Write()
hNumerOriginal.Write()
hDenomWeighted.Write()
hNumerWeighted.Write()

hAcceptanceOriginal.Write()
hAcceptanceWeighted.Write()

hAcceptanceOriginal_Transpose.Write()
hAcceptanceWeighted_Transpose.Write()


outputFile.Close()
fRef.Close()

# ============================================================
# Save acceptance results to TXT
# ============================================================

txtFileName = "Acceptance_reweighted_v2.txt"

with open(txtFileName, "w") as txt:

    txt.write("============================================================\n")
    txt.write("J/psi Acceptance Results\n")
    txt.write("============================================================\n")
    txt.write("\n")

    txt.write("Reference histogram : hPromptPtY\n")
    txt.write("Fine pT bins        : 68\n")
    txt.write("Fine |y| bins       : 12\n")
    txt.write("Analysis pT bins    : 34\n")
    txt.write("Analysis |y| bins   : 6\n")
    txt.write("\n")

    txt.write("Total events        : {}\n".format(nTotalEvents))
    txt.write("Total GEN J/psi     : {}\n".format(nGenJpsi))
    txt.write("Selected J/psi      : {}\n".format(nSelected))
    txt.write("Accepted J/psi      : {}\n".format(nAccepted))
    txt.write("Sum weights         : {:.10f}\n".format(sumWeights))
    txt.write("Accepted sum weights: {:.10f}\n".format(sumWeightsAccepted))
    txt.write("\n")

    txt.write(
        "Integrated original acceptance : {:.10f}\n".format(
            integratedAcceptanceOriginal
        )
    )

    txt.write(
        "Integrated weighted acceptance : {:.10f}\n".format(
            integratedAcceptanceWeighted
        )
    )

    txt.write(
        "Integrated difference          : {:.10f}\n".format(
            integratedAcceptanceWeighted
            -
            integratedAcceptanceOriginal
        )
    )

    txt.write("\n")
    txt.write("============================================================\n")
    txt.write("Bin-by-bin acceptance\n")
    txt.write("============================================================\n")

    txt.write(
        "{:<8} {:<8} {:>14} {:>14} {:>14} {:>14} {:>14} {:>14}\n".format(
            "pT_min",
            "pT_max",
            "y_min",
            "y_max",
            "Denom",
            "Numer",
            "Acc_orig",
            "Acc_weight"
        )
    )

    txt.write("-" * 110 + "\n")

    for ix in range(
        1,
        hAcceptanceWeighted.GetNbinsX() + 1
    ):

        ptMin = hAcceptanceWeighted.GetXaxis().GetBinLowEdge(ix)
        ptMax = hAcceptanceWeighted.GetXaxis().GetBinUpEdge(ix)

        for iy in range(
            1,
            hAcceptanceWeighted.GetNbinsY() + 1
        ):

            yMin = hAcceptanceWeighted.GetYaxis().GetBinLowEdge(iy)
            yMax = hAcceptanceWeighted.GetYaxis().GetBinUpEdge(iy)

            denomOriginal = hDenomOriginal.GetBinContent(
                ix,
                iy
            )

            numerOriginal = hNumerOriginal.GetBinContent(
                ix,
                iy
            )

            denomWeighted = hDenomWeighted.GetBinContent(
                ix,
                iy
            )

            numerWeighted = hNumerWeighted.GetBinContent(
                ix,
                iy
            )

            acceptanceOriginal = hAcceptanceOriginal.GetBinContent(
                ix,
                iy
            )

            acceptanceWeighted = hAcceptanceWeighted.GetBinContent(
                ix,
                iy
            )

            difference = (
                acceptanceWeighted
                -
                acceptanceOriginal
            )

            txt.write(
                "{:<8.2f} {:<8.2f} {:>8.2f} {:>8.2f} "
                "{:>14.6f} {:>14.6f} {:>14.8f} {:>14.8f}\n".format(
                    ptMin,
                    ptMax,
                    yMin,
                    yMax,
                    denomWeighted,
                    numerWeighted,
                    acceptanceOriginal,
                    acceptanceWeighted
                )
            )

    txt.write("\n")
    txt.write("============================================================\n")
    txt.write("Detailed bin information\n")
    txt.write("============================================================\n")

    for ix in range(
        1,
        hAcceptanceWeighted.GetNbinsX() + 1
    ):

        ptMin = hAcceptanceWeighted.GetXaxis().GetBinLowEdge(ix)
        ptMax = hAcceptanceWeighted.GetXaxis().GetBinUpEdge(ix)

        for iy in range(
            1,
            hAcceptanceWeighted.GetNbinsY() + 1
        ):

            yMin = hAcceptanceWeighted.GetYaxis().GetBinLowEdge(iy)
            yMax = hAcceptanceWeighted.GetYaxis().GetBinUpEdge(iy)

            denomOriginal = hDenomOriginal.GetBinContent(
                ix,
                iy
            )

            numerOriginal = hNumerOriginal.GetBinContent(
                ix,
                iy
            )

            denomWeighted = hDenomWeighted.GetBinContent(
                ix,
                iy
            )

            numerWeighted = hNumerWeighted.GetBinContent(
                ix,
                iy
            )

            acceptanceOriginal = hAcceptanceOriginal.GetBinContent(
                ix,
                iy
            )

            acceptanceWeighted = hAcceptanceWeighted.GetBinContent(
                ix,
                iy
            )

            difference = (
                acceptanceWeighted
                -
                acceptanceOriginal
            )

            relativeDifference = 0.0

            if acceptanceOriginal > 0:

                relativeDifference = (
                    difference
                    /
                    acceptanceOriginal
                )

            txt.write(
                "pT = [{:.2f}, {:.2f}) GeV, "
                "|y| = [{:.2f}, {:.2f})\n".format(
                    ptMin,
                    ptMax,
                    yMin,
                    yMax
                )
            )

            txt.write(
                "  Original : "
                "denom = {:.6f}, "
                "numer = {:.6f}, "
                "acceptance = {:.10f}\n".format(
                    denomOriginal,
                    numerOriginal,
                    acceptanceOriginal
                )
            )

            txt.write(
                "  Weighted : "
                "denom = {:.6f}, "
                "numer = {:.6f}, "
                "acceptance = {:.10f}\n".format(
                    denomWeighted,
                    numerWeighted,
                    acceptanceWeighted
                )
            )

            txt.write(
                "  Difference = {:.10f}\n".format(
                    difference
                )
            )

            txt.write(
                "  Relative difference = {:.6f}%\n".format(
                    relativeDifference * 100.0
                )
            )

            txt.write("\n")


print()
print("============================================================")
print("Output")
print("============================================================")
print("ROOT file:")
print("  Acceptance_reweighted_v2.root")
print()
print("TXT file:")
print("  Acceptance_reweighted_v2.txt")
print("============================================================")
