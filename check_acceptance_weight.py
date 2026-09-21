import ROOT

ROOT.gROOT.SetBatch(True)

file = ROOT.TFile.Open("Acceptance_reweighted.root")

hDenom = file.Get("hDenomWeighted")
hNumer = file.Get("hNumerWeighted")
hAcceptance = file.Get("hAcceptanceWeighted")

# ============================================================
# Bins to check
# Each entry is ONE independent bin
# ============================================================

binsToCheck = [
    (10.0, 0.0),
    (12.0, 0.0),
    (15.0, 0.6),
    (18.0, 1.2),
    (20.0, 0.9),
    (25.0, 1.8),
    (30.0, 0.0),
    (38.0, 0.6),
    (50.0, 1.2),
    (75.0, 1.8),
]

# ============================================================
# Check each bin
# ============================================================

print("==========================================================================")
print("Check weighted acceptance for multiple independent bins")
print("==========================================================================")

print(
    f"{'pT bin':>15} "
    f"{'|y| bin':>15} "
    f"{'Weighted Denom':>18} "
    f"{'Weighted Numer':>18} "
    f"{'N/D':>14} "
    f"{'Map':>14} "
    f"{'Difference':>14}"
)

print("--------------------------------------------------------------------------")

for ptMin, yMin in binsToCheck:

    binX = hAcceptance.GetXaxis().FindBin(ptMin)
    binY = hAcceptance.GetYaxis().FindBin(yMin)

    # Actual bin boundaries
    ptLow = hAcceptance.GetXaxis().GetBinLowEdge(binX)
    ptHigh = hAcceptance.GetXaxis().GetBinUpEdge(binX)

    yLow = hAcceptance.GetYaxis().GetBinLowEdge(binY)
    yHigh = hAcceptance.GetYaxis().GetBinUpEdge(binY)

    # Weighted numerator and denominator
    denom = hDenom.GetBinContent(binX, binY)
    numer = hNumer.GetBinContent(binX, binY)

    # Calculate acceptance directly from weighted histograms
    if denom > 0:
        acceptanceFromWeightedHist = numer / denom
    else:
        acceptanceFromWeightedHist = 0.0

    # Acceptance stored in acceptance histogram
    acceptanceMap = hAcceptance.GetBinContent(binX, binY)

    # Difference
    difference = acceptanceFromWeightedHist - acceptanceMap

    print(
        f"{ptLow:.1f}-{ptHigh:.1f}".rjust(15),
        f"{yLow:.1f}-{yHigh:.1f}".rjust(15),
        f"{denom:18.6f}",
        f"{numer:18.6f}",
        f"{acceptanceFromWeightedHist:14.8f}",
        f"{acceptanceMap:14.8f}",
        f"{difference:14.3e}"
    )

print("==========================================================================")

file.Close()
