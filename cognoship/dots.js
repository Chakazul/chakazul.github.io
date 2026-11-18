"use strict";
// ============================================================================================
//  The dots channel of the game ([dot:<n>|off] in script.txt), ported from Carlman: a separate
//  Lenia field of small stationary solitons (Carlman's dot: rule0, mu 0.38, sigma 0.07, at half
//  the kernel radius) that breathe under their own rule but never move. Any creature eats a dot
//  by covering it (its mass erases the dot's, see shaders/dotsim.glsl); a dot bitten deep enough
//  dissolves, a graze heals. The game keeps <n> dots on the board: an eaten one is replaced at
//  once, at a random spot away from the creatures and the other dots. The player's ship eating
//  one scores 10 (the score, shown at the top, starts from 0 when the dots are turned on).
// ============================================================================================

const Dots = (function () {
  const C = Core;
  // Carlman's dot soliton (models/solitons_direction.json rule0_mu0.3800_s0.0700_R18), 22 x 22 at R 18.
  const SOLITON = [
    [0,0,0,0,0,0,0,0.0015,0.0208,0.0451,0.0599,0.0574,0.0391,0.0142,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0.0068,0.0697,0.1637,0.2481,0.6059,0.6721,0.6601,0.5712,0.2241,0.1349,0.0475,0.0009,0,0,0,0,0],
    [0,0,0,0,0.0384,0.1698,0.64,0.9472,1,1,1,1,1,1,0.8626,0.5417,0.1276,0.0184,0,0,0,0],
    [0,0,0,0.0551,0.2293,0.8843,1,1,1,1,1,1,1,1,1,1,0.7632,0.178,0.0268,0,0,0],
    [0,0,0.0384,0.2293,0.9695,1,1,1,1,1,1,1,1,1,1,1,1,0.8295,0.1731,0.0143,0,0],
    [0,0.0068,0.1698,0.8843,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.7354,0.1142,0,0],
    [0,0.0697,0.64,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.3662,0.0333,0],
    [0.0015,0.1637,0.9472,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.7858,0.1061,0],
    [0.0208,0.2481,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.1836,0.0035],
    [0.0451,0.6059,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.2408,0.0167],
    [0.0599,0.6721,1,1,1,1,1,1,1,1,0.9983,0.9991,1,1,1,1,1,1,1,1,0.5115,0.0258],
    [0.0574,0.6601,1,1,1,1,1,1,1,1,0.9991,0.9999,1,1,1,1,1,1,1,1,0.501,0.0245],
    [0.0391,0.5712,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.2268,0.013],
    [0.0142,0.2241,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.9638,0.1633,0.0014],
    [0,0.1349,0.8626,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.7105,0.0841,0],
    [0,0.0475,0.5417,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.2158,0.0201,0],
    [0,0.0009,0.1276,0.7632,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.6274,0.0814,0,0],
    [0,0,0.0184,0.178,0.8295,1,1,1,1,1,1,1,1,1,1,1,1,0.7086,0.1279,0.0038,0,0],
    [0,0,0,0.0268,0.1731,0.7354,1,1,1,1,1,1,1,1,1,1,0.6274,0.1279,0.0104,0,0,0],
    [0,0,0,0,0.0143,0.1142,0.3662,0.7858,1,1,1,1,1,0.9638,0.7105,0.2158,0.0814,0.0038,0,0,0,0],
    [0,0,0,0,0,0,0.0333,0.1061,0.1836,0.2408,0.5115,0.501,0.2268,0.1633,0.0841,0.0201,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0.0035,0.0167,0.0258,0.0245,0.013,0.0014,0,0,0,0,0,0,0,0]
  ];
  const SOLITON_R = 18;
  const RULE = { mu: 0.38, sigma: 0.07, dt: 0.1 };   // Carlman's dt, per step
  // Patterns that [dot:<name>:<how>] puts in the dots channel instead of the round dots, one at a
  // time and not replaced when eaten (a beacon). The channel runs one rule at a time: the dots' or
  // the pattern's. From Lenia/Python/animals.json (Lenia.py's numbering: kn 1 = exp core, gn 1 =
  // Gaussian growth, as dotsim.glsl's).
  //   Heptadentium scintillans (2D7s): checked in numpy with the dots' maths (exp core, rings
  //   1,5/6, Gaussian growth, dt 0.1): at R 18 it stays put (0.4 cells in 3,000 steps) and its mass
  //   holds within 1 %; at R 12 too; at R 9 (the round dots' R) it dies.
  //   The others below run as animals.json means them (LeniaND.py: kn 1 = the quad4 core, gn 1 = the
  //   polynomial growth, growth: 'poly'), each checked in numpy with those maths for 3,000 steps
  //   (its mass within 0.67-1.05 of its own; Tetrahelicium at R 12 too). Most stay put; Circoechinium
  //   and Gyrogeminium wander (about 25 cells in 3,000 steps): a beacon's site follows its pattern (frame).
  const PATTERNS = {
    heptadentium: { name: 'Heptadentium scintillans', R: 18, betas: [1, 5 / 6], core: 'exp', mu: 0.22, sigma: 0.025, dt: 0.1,
      rle: '22.NpFK$11.wQwSuLrA6.pLtLwX2yOyGpB$9.tE3yOyIvOrN5.rPvMxT3yOyLuV$8.tH2yOtAqUtIwIwEtIpN2.MtJwGwXxV2yOyLuSpW$7.qCwNwWrA3.pNtPvTyFxHxN3yOyJ3yOyNuMpC$7.qWvHvArI5.pXtOxM9yOvXpC$7.pKtOvTvPvHtGrC4.pAtL8yOyMrK$7.AqS3yOyHwWtV4.XwG3yOyIvSsUyCyOwN$8.6yOyCsS3.BuP3yOxDrN2.wQyOrO$6.AvU7yOvJ4.tN2yOyLuM3.pEyOyH3.E$5.qQvGxRvLtVvJyI3yOvI4.rA2yOxDrB3.NxAyOtLSqHrQrAA$2.pVqDsLwWyOsBODSrSxFyOyN6.qTxRuD4.qBwV2yOyJuVwFwNrT$.uGtWuTxMyOuE5.2I13.XtUxP4yOyH2yO$pQxGxKyH2yON19.OtAxWyM3yOyJvR2yOpX$vD4yOxQqM19.OxD4yOyCtVpTwEyOxA$xU4yOxCvDtTsJqL17.wR3yOxStD2.sHyOxI$xJyOxLvFxRyCyK2yOyHxC17.xFxNwIsX2.DsVyNvB$uPyOvTuTxT6yOwH16.qIsAqV3.qTvIwHrU$rGwQvTvVyN7yO21.pNvTtFrXD$MtBuOvK4yOyEyHyArX20.pAxAuPW$.pAqRrR3yOuLrUqBM20.pHxDyOA$4.wCyOwHpH19.pPqKqPrQwEyOtG$4.pU2yO19.sRwUwNwAxJ2yO$4.MyFyO18.sJxJ6yOP$4.pRxTyO3.VrTrR12.qTuUyE5yOtCqHE$4.rVyDyO2.qXvCyIyOyApW11.rTvGxQyKyMyIxLvNuFqL$3.qFuLyOyApAqXvE4yOsH6.rIpK3.HrVuLvQvIvAtMvFwPsQ$3.tExAyOuDqLtFxJ3yOwVqT3.pGwNyOyAsUJ3.pBsRuMtUtLsWwWyOsJ$3.vMyJyOsXpGuFyB3yOvIU3.rQxR2yOvWqX3.qAxE2yOyA3yOK$4.2yOtQpIsFyA3yOuIA3.sDxS2yOyIuUqWpBqU7yOtP$4.vWyOxCrRsNvGuU2yOvWB3.rXxL8yOyKxJxMyG2yO$5.xKyOvOuEuGsEpF2yOqJ3.rUxE6yOyDtDqIqRsRvAwM$6.tJvIuCsDpJ2.yJyOqB2.tBxP4yOxVpO5.qV$7.PpKP4.vIxUsWtAyL3yOwIuBqOB$15.uFyF4yOxWvVuCqTB$15.rMwFyM4yOwLuAqH$16.sWwV4yOwFrXJ$18.sBvOvMtCpR!' },
    'helicium': { name: 'Helicium serratus turbulentus', code: '4H3t', R: 18, betas: [1 / 4, 1 / 2, 3 / 4, 1], core: 'quad4', growth: 'poly', mu: 0.2, sigma: 0.02, dt: 0.1,
      rle: '53.J$42.pN8.GLG$42.qJ8.pIG.T$13.2LG26.E3.sC4.G2.pN$12.JQOLQpDpIpDL20.V3.LsP7.pPpUpXV$11.LpDVOVpN2pXpPpAJ19.OQ2.qJpU.2BO3.pSqWrLrG$9.LpUqMqJqCqEqMqTqOqEpP2pNpFTLE15.EQpPrGrQ2pDQsKsF3.qCsCrLrS$7.EpKqMqOqJqMqRrArGrNrIqJqEqMqOqMpUpAJB12.ELQpDqHrIrVpUpIpXvMsHBE.qJsFqCrSV$6.EpAqH2rGrDrIrNrQrGrNrVrNrLrNrLrDqJpNVLG7.E2QpFpNpPpXqEqJqOrG2rVtQuCqC2pXpUpNpAOsCrQO$6.BVqMrSrVrSsCsWtGsUsHsKsPtBsUsCrIqWqMqCpNpFpDrGsHtOsW2sHrSrNrDqWqO2qM2qHqErDsRuNtJqCpNpPqHqOqCpKpDsPtJsHqH$6.JqErVrQrIsCtBtVuCtOtBtGtQtVtQsWsC2rVrStJvJxMyIxRwQvEuFtJsPrXrIqTqJqEpNpIpPpNqEqOqRqEpApNrQxRtGpUpNpItJuCuFsW$5.OqOrLpXTpIqRrLrXtGuKuIuPuXuNuItQ3tGxC3yOwQtTrVqHpNpDpA2pDpAVOGJOJE.OpXqMrSwSyLqOqC2pXuNuXuPtV$4.pAqWrApULOqErVsCrXtDuXvTwIwDvJuSuNuPuU3yOyDsPpF17.JqHsMtJsCpDEGTpXuI2uCuI$3.TqRrLrSrIqOrLuAvRwA2vOwDwQwXwSwDvOvRwD3yOxMqM19.OpPqHpP4.BLpKqEtOvBpS$2.JpXrArS2sRrNsMvBwQ2xEwVxE2xPxE2wQxC4yOqJ18.BJTpPqEqC2.qCtBLJLpKuAwDtTqW$.BpApXqTrVsMqHLpKrQsPuAwSyByIyOyGxHwVxJ4yOsFE17.EVpFpKqEqTrAqOrX2wVJ.EpKuP2wAtT$.EVqErNrSpA3.qEqJrQvW4yOxRxEyI3yOxMqE15.pA.EpDpSqEqJqWrSsFsPrXsMpK3.pPuPvBvTuN$2.VrIrSpA3.rX2vTwDyL4yOyLyB4yOuST15.qTpUpIpUqOrGrLrVsPtBsUV5.pDsWtJvOtQ$2.qCsMrApKG.qTwX14yOtO16.sCsKqRqWrNsFsKsRtJtTtGQ.2L4.qJvWsK$.GrXsKrVrXqTOqEuPvJuPwNyL2yOyL7yOrG16.tDuSsMsCsMtDtJtLuFuPuK2sR2wIG3.pUvOsWpI$.qEsPsMsKrXQ3.QpNtBxRyDyBxUxW2yDyByI2yO17.uFwIuItLtQuFuIuFuXvOvR2vGxHuP4.rNuNuFrXO$BrGsWtBsHpI4.qTrVuPxMxRxEwVwXxEwXwVxMyIxR17.vMxJvTuPuS2vGuXvO2wIuIrDqH4.QsKtGuAsM$LqJsPtVrG4.uUxUwVwSxExCwIvWwAwDvTvRwNxWwQ17.wLyBwVvWwAwIwFvWwIwXwQsUJ5.rAsFsMtOrD$.LsHuNpN3.qCxCxHuSvE2vWvMuX2vBuPuKvGxEvM17.xMyIxUxCxExMxEwXxJxUxMvJsKtBtTpP.BqJrAsHtJpA$2.tJuCqRpS3.pPpIEqOuUvBuP3uAtOtGtOvOuI17.yL2yOyDyGyLyBxUyB2yLyBxRyOyGrL3.pPsWrVE$2.uSuAtBsRE5.qCuIuNuAtD2tBsPrXrVtOtD16.pU7yO2yL2yOwXvTwFtL4.rDsRqTT$.TvEuUuFsR4.qRsFsKuAtVtJsM2sCrNqWqJrIsC16.sH5yOyL4yOyDtDqJpX4.qTsFrXqWpK$.pKuNvWvGrG3.qHxRwDrXsMsRsMrSrGrDqOpXpFpDqOO14.TtG4yOxWxR5yOuSsHrQpP2.rDsRsCrSqWQ$.JrIuXwIpP3.qRuApU.L2rQrAqJqEpSpDE.LJ14.qEvW4yOxPxU5yOxWxMwVuF2qW2sPsCrSpX$3.sRwNrDpAE6.rArDqMpPpFpAE17.JsC4yOxJxCxPyLyOyIxUxJwXxEwLtOqHpNqOrLrVqWE$3.sHwDuKuUqRE3.EqErXrApXpAJB18.qMxR3yOxCwIwFwXxRxWxMwIuF2sPrQpA2.pNrLrGQ$3.sPvWvMvGrGpKQVuCxPuIrIqEpFL19.qTxC3yOwAvOvJvTwL2wSwLvWtOrGqRpSLEpUrSrAV$3.rSvEvGuSqHpXpUqHxUwIQ.BTG17.pKuC4yOuPuKuFuNvGvTvWvJuUvGvEuItLsFrGrSsMrDpFG$3.pKsPuFuNqEqCqE2qHOB.GrAqTpSJL2JQTpDpKpPpNpIpNpXqTrVtQxC3yOwItB2tDtGtTuNvBuSuFuAuIuSuNtLsPsKsPrSqEpAG$5.rVtQpPqCqJqEpKVpFqMtQuSsWqTpApFpIpSqEqCqJrArLsCsPtJuFvJwFxHxWwXvJrX2rLrQsCsPtDtOtTtOtBsFrXsMsRsC2rIrQrIqHV$5.pItBrAtBrQpNpPpUsFxEuFqOrGrAqMqOqJ2qO3qMqOqTrGrLsCsRsCrSrDVpDpNpXqEqOrIsF2sK2sFsMrXrGrDqWqMpXqHqWqMpD$6.sWsKtBrD.GJuCwFpU2qHsPsMrNqM2pUpPpFQJB8.EJVpNqCqTrNrQrLrDrGrXsKsFrNrD2qTrApXJ$6.sKsHsCqH2.GrNpPpU2qCtDsMqTpSpFOG14.JpFpPpXqHqMqOqRrDrLrIqO2qEqMqET$6.qHqWrIpU2.G2TQGsRrG.OpXpKL15.BOpAVQVpNqMqOqRqMpUpN2pSV$8.pPpU2.Q3.pPtV3.pUpS18.BEBGVpPqCpSpPpNpPpSpNL$9.pU.qJpS3.VrI3.pX22.G2LJGJQpDVJ$9.pNpFpUpD8.qM3.E$10.pPVL7.pAqH$25.L!' },
    'tetrahelicium': { name: 'Tetrahelicium serratus turbulentus', code: '2H4t', R: 10, betas: [1 / 2, 1], core: 'quad4', growth: 'poly', mu: 0.29, sigma: 0.043, dt: 0.1,
      rle: '38.pKVLO$38.qRpS2pFpA$38.rDrAqWqTqH$37.VtJuC2tQtBpX$37.qWuFrQrIwFvWvEpS$8.OTJ24.LpPrIsCsHsCxUwVwFwDsH$6.JpPqRqWpXQJQ2pFTOJ14.pFpX2qJpXrArVuFyBxMxCwQwI$5.BpNrDsCrXrIrArLsCsMsPsHrQqTpKB10.pFrArQrApPJGpSuU4yOxWwX$5.pKrVtOtQsMsFtQuXvT2wDvWvEtVsMqJJ2.E2qJpFpApKqCrLsUsRqHE3.rIsUpNqRwLyOyI$3.EqWtOuKvTwLtOtVwLxP3yOwNuXvGvMtTrSsRuNwFwXvRuIsUrQrDqMpXqCT4.B2.EwD2yOqM$2.VrQtLsPrLuKyOxCyB4yOwXsMrIuNyD5yOyLuPrApAG8.GpD2.rS2yOyLvGL$.VrGtBtTsMrQvW7yOsK2.tDxP6yOrI9.BpDpXqWsFuAyD3yOxJqT$.qCrStLvOwDwQ7yOxP3.rAuIwAxE4yOsH7.LpUqErDrXtBuK2vGwIwXxHxErA$BqCrGtJxP8yOyIsC4.TBrD4yOvJ7.qTsMsPtGtVuUvRqT.pPtDwLvWpF$.pUsUxEyLuAuXyLyIyBxR2xJwD7.rA2yOxRyOuA7.tLvBvEvMvWwNwVE.qEvTwQsF$.sFwQvTqE.EwVwNvWvMvEvBtL7.uAyOxR2yOrA7.wD2xJxRyByIyLuXuAyLxEsUpU$pFvWwLtDpP.qTvRuUtVtGsPsMqT7.vJ4yOrDBT4.sCyI8yOxPtJrGqCB$rAxExHwXwI2vGuKtBrXrDqEpUL7.sH4yOxEwAuIrA3.xP7yOwQwDvOtLrSqC$qTxJ3yOyDuAsFqWpXpDB9.rI6yOxPtD2.sK7yOvWrQsMtTtBrGV$LvGyL2yOrS2.pDG8.GpArAuPyL5yOyDuNrIsMwX4yOyBxCyOuKrLsPtLrQV$.qM2yOwDE2.B4.TqCpXqMrDrQsUuIvRwXwFuNsRrStTvMvGuXwN3yOxPwLtVtOwLvTuKtOqWE$2.yIyOwLqRpNsUrI3.EqHsRsUrLqCpKpApF2qJE2.JqJsMtVvEvW2wDvTuXtQsFsMtQtOrVpK$2.wXxW4yOuUpSGJpPrArQrApF10.BpKqTrQsHsPsMsCrLrArIrXsCrDpNB$2.wIwQxCxMyBuFrVrApX2qJpXpF14.JOT2pFQJQpXqWqRpPJ$2.sHwDwFwVxUsCsHsCrIpPL24.JTO$3.pSvEvWwFrIrQuFqW$4.pXtB2tQuCtJV$5.qHqTqWrArD$5.pA2pFpSqR$6.OLVpK!' },
    'circoechinium': { name: 'Circoechinium ventilans', code: '3ECv', R: 18, betas: [1, 1, 1], core: 'quad4', growth: 'poly', mu: 0.29, sigma: 0.0345, dt: 0.1,
      rle: '24.pNuKwVvGpP$24.qWvJ3yOyIsF3.qRtVtD$15.pSuFwSxPyD3yOxE2vEyD5yOwLrNsMwFxUuK$14.QuF18yOvMwS2yOuK$14.qHuKyG22yOxC$13.sCwQ26yOrApFrDqE$11.tV10yOyLyGxMwS16yOuAuUvJrI$8.pSvE12yO2yGuX2.tGyI2yOyLvOwI8yOyIwLxEuN$7.rIwX9yOyDyOyDxWyIyOyLV2.vM2yOxMqH.sM10yOxCqC$6.pDtD8yOvTpI3.qRxC2yOtL3.vOyOvJ2.G11yOuS$6.VsM7yOyDwSuFqE5.wLuI4.GqH3.3yOwXxM6yOyLrA$6.pPuX7yOxExJ2yOsP5.pS9.sRyOvWqWpNuI6yOxUQ$5.pFwI7yOvJsWrGxJ2yO16.BL2.pN7yOwFqTB$4.V6yOwVuFqM4.rNxR20.rG8yOvTuAO$4.6yOwIqRqWV25.GuN8yOxEwDsF$3.wA6yOxEuFxRyOqC25.sCyOyIxC5yOyIvTsK$2.qJ7yOwXvJ3yOE25.VrNrDtT5yOwQrI$2.rI7yOuS.EwQyOwF28.OwF5yOsF$2.pKwA5yOyBV33.pIuU5yOvT$2.EvG5yOtT17.GB14.rDuX6yOxWT$2.pK4yOwVvEtQO13.rLxJ4yOtT11.EuS7yOyItD$.VvO4yOqWsRyDyOvE11.xH7yOxUpS10.BvByOyL4yOyLxCJ$.uU5yOtOsR3yOxR9.11yOT11.vOwX6yOsR$.6yOwXrIqOtGuKrL8.wA3yO4.vB3yOwI11.pIuF5yOyGuNB$O7yOpU12.3yO6.uF3yOqH11.rAxU4yOvWtTQ$pN7yOpK11.rX2yOtJ7.3yOuA9.qHrSrVvJ4yOtOrD$E7yOqT11.tG2yOpD7.xH2yOvW8.VvTwSvEwD4yOsR$.rA5yOvJ2tLqT9.sR2yOqJ7.wV2yOvR9.vByOxRxU4yOuA$2.uS3yOsWrAvG3yOpD7.qO2yOwN7.3yOtJ10.wLyByD4yOwI$2.rS3yOqRpFuC2yOuI9.3yOG5.rL3yOpF10.sRyG5yOwX$2.tO4yOuAsKqT11.vM3yOrG3.sC3yOvO11.sC6yOxM$.pDyL5yOvM13.yD9yOxW11.qMtL6yOyIQ$.sR7yO14.wL7yOwA10.rNvTuXtJxH6yOpI$.tB7yOrX14.qOvR3yOwLqO11.rGxHvEsHuK5yOvTO$.pK7yOsUuP2yOxJ29.wAuPsHuK4yOtDrS$2.yI5yOuSqWvO2yOpD29.vWwIvGxU3yOsMBE$2.TvM4yOrApNtVxRsC29.B6yOyLQ$4.qM3yOxMvJvWuP27.rQpA.qO6yOyI$5.uP6yO4.tVqH21.wVvRrIsP6yOyD$5.qE6yOwQ.rN2yOtQ21.wNuUqRqT6yOyB$5.rG7yOsFwD2yO22.wAtLpPpS6yOxH$5.sM7yOrNuFyOwL4.tDuX10.rIuN3.rQyLvOsHtL6yOsF$5.pK7yOqOsPvWtJ3.qR2yO4.sUtT4.uXyOrD2.7yO2yLyOwF$6.vB6yOxE3yOtL2.3yO4.2yO4.wQyGsF.sM7yOuNpUrLpS$7.tT10yOvJtByByOuN3.pD2yO3.pIxRwFqO.xJ6yOyG$10.rQ8yOuPtVvEsH3.vJyOxE2.qJyGyOwQrIrG7yOsM$11.qH7yOxRtJvGwQwNwVxCxMxCsWpPtT4yOwX8yO$12.wN14yOxHtQwF13yOtO$12.tG25yOyIyL2yOwN$13.wS23yOyLsRuCwXuS$14.rDvTwVvWvOxW16yOwQL.pKqE$19.pFtL14yOvG$20.pUvJ5yOtVrGrVvM2yOyGqT$21.GqErDqRT3.LrDsKqC!' },
    'gyrogeminium': { name: 'Gyrogeminium serratus', code: '3GG2s', R: 18, betas: [1 / 2, 1, 1 / 2], core: 'quad4', growth: 'poly', mu: 0.25, sigma: 0.031, dt: 0.1,
      rle: '20.AHMRSROLIFB$12.BF2GFEFJWpNqCqIqHqCpSpIXQLGEC$11.GNXpHpOpWqCqDqHqPrCrLrPrMrAqKpVpMpGpBTPMKHDA$11.EMpCpXqSrIrUsFsKsMsLsJsErUrHqT2qNqMqIpXpMpGpBVQKEA$10.ABGWqGrOsIsPtCtRtVtIsMrWrQ2rPrUsEsKsIrUrDqKqApUpTpQpKXMD$7.MpMpXqApVpWqLrNsStOtWuGuRuWuFtDsNsPtAtItPtWuEuHuEtQsXsLsJsOsSsRsErIqJpNUG$6.pOrHsEsOsXtItWuKuUvBvJvPvS2vUvOvHvGvM2vR2vOvSwDwOwWwXwTwOwMwKwCvLuNtIsDrAqAXE$5.pLrMsMsUsTsOsKsMtLvEwTxHxDwW2wRxDxSyIyMyGxWxRxQxWyH6yOyExPxCwLvMuHtBrUqPpIH$4.GqArLsIsVtItGsLrNrQtOwQ2yOxXxQxW9yOyJxWxJxGxQyGyJyFxWxUxTxJwJvCtTsLrApKG$4.FpHpXqKrGsQuDuKtUtRuTwV8yOyNyAxTxSxNxBwDvAuAtVuSwGxGxSyCyL3yOxOwFuTtFrIpLF$4.pLqXrTsLtKuUwEwSwWxCxMyC6yOyNyJxSwOvFtRsPsArQrMrPsGtLuMvMwGwVxOyH3yOyBwSvKtOrKpKE$3.qQtUvPwLwTwUxAxN4yOyHyEyFyJyMyKyGxWxKwDuArKpE2.CpKqNrArPsMtMuIvCvVwHwXxUyDxRxCwJvFtIrDpFC$2.pWuMvQuStErDpFDpEtHyEyGyCxR2xMxRxVxPwRvCtWsEpQ6.QpFpQqJrHsGtDtOtUuPvQwJwMwIwEvPuIsJqJS$2.tAwAvLuHsPqE3.XuG2wXwSwRwUxCxDvLsBpQB11.FpFqKrQsVuFvDwEwTwUwOwGvTuQtArCpJE$.pAuBwFwHwCvSuLrSpNpGrBtXvOvPvOvSwEwNwJtEU15.pPrCsOuBvQxEyEyNxWxCwEuTtErJpQJ$.UsOuPvEvNwRyAxKvRuRvDuWuSuOuLuPvEvMuVrB16.pHqTsHtSvKxCyH2yOyDwPuVtBrHpRL$.RsPwByK4yOyIwQvIuOuEtTtMtLtQtMqR17.RqErOtCuRwNyB3yOxLvKtErDpMJ$.rSxG2yOwNsOqDqGsVuNtUtLtEtAsTsLrWM19.pDqNsCtUvQxLyL2yOyCwHtTrIpJF$pJwF2yOuMH4.NrUrSrTsBsFrRpG21.pQrMtMvMxEyDyNyMyDwWuQrXpLD$sI3yOwFqWC4.pPqHqKqXrNrD22.qJtAvMxNyM3yOyGxFvHsOpTD$uG5yOwNtIrDqEqGpQpIpHpUqNqH21.pQtDwGyH6yOxJvMsWqBD$sKwSyByJ5yOuCqRpNTPXpGR6.BGIF10.UtOxCyN7yOwWvFsSqAD$LtSxX6yOtRqDTFHK4.DILSpGpRpWpJJ8.JvCyI7yOyCwLvDtVrXpOB$DvLxW2yOtIrDqKpGBA2.IE3.JpGpTqBqIqUrCrDqRpIG2.CEA.BwU8yOyHvVtXsNrNqLW$rGvNwUwDsFrHqPpWS3.HpIpBGBKpQqTrDrErHrDqVrCrJqPpPpJpTqKqQqDpD10yOxIuIsAqGpHSE$tNuVvQvMuOsNrUrCpTPINqFsBrPpKpEpVrDrKrCqXqNpVpKpWqUrFrErJrSsBsOsW11yOwVtDqJRC$uFuLuTvMwLxFxPyExMtGrBrDsWuJtVrMqOrDrHqMpVpMXC.GpGpXqKqOqKqIrLtO11yOwFrWpB$rOtGtWuNuWvCvDvRxCuFsErTsTtTtSsFqXrAqHpCI6.BIH2.pAvFyEyNyMxTwOvUvNwDyL2yOuRqJD$.UtVvBvLvIvCtOrMqWqKqXrHrVsAqQqJqEpC14.uFwAwUwGtNqL3.yAyOwVsAQ$.qDuBvNuJqDpHXpJpWqIqUrDqRqDpPpHV15.rBsOtDqU6.uNrWU$.qWsOrCpPqA3pXqAqEqHqApIRH18.PH7.PB$BqApUpDpLqCqVrJrFqDpK2pIQ$pFpQpXqOrDrMrSsBsPpE$MpOpKpLpNpOpRqEQ$2.RpOpWpXJ$2.UpXqC$3.G$7.HP!' },
    'decadentium': { name: 'Decadentium rotans bioculi', code: '2D10rb', R: 18, betas: [1, 1], core: 'quad4', growth: 'poly', mu: 0.28, sigma: 0.036, dt: 0.1,
      rle: '22.ETJ$21.pXtDvJvMtBpI$19.OtD3yOyBuSpI6.JpIpUpXpNpAB$18.JtL4yOyLvJpU4.GrDuPxR2yOxWuSrG$18.qWxH5yOxPtG4.qJvM5yOyIwNsF$10.J7.rVxU6yOyBuPqC2.rSxR5yOxHtLpP$8.LrG2tLrV5.qMwQyOyIxPyG4yOxUwIvW7yOwNqW$8.qRvWyGwLsUB4.LuA2yOxRvRuFuKyG11yOxHqO$7.LuK2yOxRvRtBqC2.pAuF3yOyLvTrD2.pPvR5yOyGxH3yOtB$7.qM5yOyIxU2xJyD5yOwXqT3.pXxC4yOwQrXpPyDyOxHJ$7.tD15yOyDrQ4.tV4yOuA3.wNyOuX2.pDrArNqCE$6.EuA4yOxEtVpNBLpUsMxC4yOtG4.rL3yOyBrG3.qT2yOuApKsUwDyBxPtVpP$6.BsU2yOyLwIsK6.qHwI3yOvG4.pP3yOwQT3.pX3yOxRwS4yOwFqJ$7.qEuU2yOxUuSqJ6.pUwI2yOwV4.BxC2yOuU4.rX10yOvOpU$8.pXuA2yOyIwLtDpN5.LtGyLwI5.pSwXwSrD3.pPwL10yOyLsH$9.pA5yOyIwIrQ5.pFrG11.TvR5yOxPyG3yOxJwFtG$9.uA7yOyGrG17.sF4yOyLvErLsF2yOxUtQrLsK$8.rL9yOtG17.sH4yOuN3.vWyOvRpI$.pXrXrA2.BtTyLyOuPtLvT5yO19.3yOvR4.qHyOtD$.sMwNuAsRtVwXyOyDpI2.QrS2yOwX20.uAxUwIpX4.sKyOqJ$JuKxWxExMyL2yO6.QB22.pXqM5.wNyL$OvT5yO36.JsHyOyL$GvR4yOvR5.qHrD8.pPV15.ErNtQuUwN3yO$.tV4yOwDsUrNsRvMxRyLyDtL6.3yO6.3yO6.tJ7yOrArIrLqEB$.qT4yOyIxWyD6yO6.3yO6.3yO6.8yOxJxMxPvWqW$.JuU12yO6.3yO6.3yO6.12yOvRO$2.pPtTvGvEvB4yOyL2yOvT6.3yO6.3yO6.xM3yOyLxEwDwNyD4yOqM$3.OpDVT4yOyLyDvGpD6.tV2yO6.rLyOB6.sFwAvWtVrGpDpArXxP4yOsF$7.2yOuKqTpNTB26.V5.pS5yOsW$6.JyOwSO29.G6.wX2yOyLyIyOsM$6.qWyOrX5.sPtG21.VyOxEQ4.qWyOyLwLvEwAwLrN$6.tOyOpK4.tDyDyOxM19.tQ3yOwIrNpFpNuIyOxErXpDpNsWsFqJ$5.QwDyOL3.rSyD3yOpA17.J6yOxE2yOwA4.QTG$3.qHpXtBxWyOvJ2.rVxP4yOrS17.tJ8yOxU$2.pIvEuSwQ3yOuXvGxW4yOxRpX11.G5.pNwS7yOB$2.pPwDxUyL9yOxPrN3.B8.uAtQQ5.sRvWxU3yOwSqTE$2.GtL11yOsP4.uNyIyOrL5.yByOxUrD6.rGuPxJ2yOvTrLJ$3.pXvO9yOxJpF3.OwV2yOxJE4.xC2yOyIrN6.qRuSxM2yOvMqH$4.pSuU3yOyGuNuA2yOxUQ3.rDyD3yOpN4.vJ3yOxMrVJ4.rSvRyB2yOyDrV$5.LrNuSvJtLqJErLyIyOqC3.tQ4yOrG4.tO5yOvTtLtVyL7yOsC$7.GpAT3.tQ2yO2.pPwD4yOtO4.sHyD15yOqT$13.BwV2yOtLuKxU4yOxEqC3.rQxE4yOyGvTuIuFvEwQyB3yOvGT$14.tL3yOyDyL5yOwVqJ.QsHwIyL3yOqC4.rLvBxRyOxJrD$14.sCxR13yOwIwVyD2yOuIV4.GtGwXvMrSO$13.LtDxM6yOuNsHtLwIyL7yOwNqM6.pKpSG$12.pSsRvTyD5yOwFqJ3.sWxH6yOwVqR$13.rGxM5yOxWtGV4.tBxJ5yOuNpA$14.qHsRuSvTvMtVrIT5.rAwAyL3yOwAqE$15.GJGB8.rLwDyI2yOuFpP$28.qEsHrXqHG$30.E!' },
    'decadentium-3': { name: 'Decadentium rotans trioculi', code: '2D10rt', R: 18, betas: [1, 5 / 6], core: 'quad4', growth: 'poly', mu: 0.28, sigma: 0.035, dt: 0.1,
      rle: '18.pDtV2yOxHsF$17.pDuP4yOvWpK7.EQG$17.qMwS5yOvOQ5.pDsMvTwQtG$17.pXwF6yOyBsK3.TtD4yOxR$9.uUqO6.pAvB3yOxCuIvE3yOxPqTsUyL4yOxRsM$7.OyLxJuFpN5.wQ3yOxWsP3.sH9yOyBrX$6.G3yOxJvBsWrVtBwN6yOuX4.rI4yOyGxC3yOsK$6.tO15yOyGsK3.pK4yOwXqOJ2yOuF$5.LwX7yOrQpKpFrLxM4yOwVpP2.qW4yOwS2.GyOxH$5.pDwN5yOxM5.qR4yOyLtV2.uC4yOxE3.vTyOpI$5.LtG4yOyGsU6.vO4yOuXVOwF4yOxC3.rG2yO2.pDrGrXpX$6.pIsRwX2yOyLuSO5.rI3yOwVO2.tG4yOvT3.pN2yOwSrGuSyB2yOtQ$7.JuK3yOyBuSrA4.LrArL6.sP2yOxWsR3.qM3yOxP5yOtJ$8.6yOyLxWwXuP12.tBrN3.OtV10yOG$7.tT10yO16.LtB11yOtT$6.tD2yOuXrIrDuK5yO16.tB4yOyDwXyB3yOyIxRxC$.pSpIpKrIvMyOyLrD3.QtJ3yO16.V4yOwXsFqWvB2yOyDvJtVtQ$.uKvBwFyD2yOT6.2pF8.3yO6.V3yOxHT2.pA2yOvBpN$TyByL3yO18.4yO6.2yOyIqH4.yOyLqE$rV4yOrQ3.pKpXJ12.4yO6.xMyDuS5.yOxW$sK3yOyGsFqJsUwLyDyLxUtJ12.2yO8.uUqH5.yOxW$qO3yOyGwFwNyI4yOyL22.qOL4.pS2yO$GtV11yO23.sWQB.pAuX2yOB$.QrNuX8yOvJ5.3yO14.EyOxJvOuUwS3yOsWsKrXpU$3.JwQ7yOL5.3yO15.10yOyIvGpD$4.qR2yOuI2sMsHpK6.3yO6.2yOyG6.12yOsH$4.pNyOwIQ10.xR2yO5.4yO6.vWyL2yOyGxPxMyG4yOvO$4.rGyOrS19.4yO7.sWtVtLsPrXtGxE4yOwX$4.uAyOpA3.sFwAvR13.3yO13.Q5yOvM$4.wSyO3.sWxW3yO28.pP3yOyDyBqM$3.sKyIyOrX.rXxU4yOwV19.tLvGrD5.2yOwSuAtOuF$2.sCwI2yOwQtJwV5yOyD18.4yOtJpDBVvRyOvBT$2.wAyI3yOxR6yOsC17.vG5yOvGtOyByOwI$2.yD9yOxWrI3.rVvGG12.xC9yO$2.rS8yOxCqH3.BwF2yOxEpK.qWuPxPrGE4.wA8yOqM$3.vR3yOxW3yOtJ4.pDxH3yOrDqOtQxH3yOvTqC2.qJuXxCyByIyL3yOpI$3.OtBuKtBqTpN2yOvO4.qRxW2yOuSO2.uAyL3yOxJqO4.TrDuPxW2yOsKL$5.LG3.2yOpD3.tByI2yOqR4.vO4yOuUT5.qOxR2yOwFqM$10.rN2yO2.pFvT3yOpK4.pDxW4yOsPL4.vB4yOsF$11.wF2yOuCvByB3yOrN5.uN5yOxUvMuAtQ6yOrG$11.tD3yOyL5yOqT4.rGxU7yOyLyIyL3yOvE$11.tGyL9yOwXV2.rGwS4yOsUG.OqErXuAwAuP$10.pKvJ6yOsHJtBxW4yOwVxW2yOwApD6.B$10.pNxM5yOxMqW3.qEwX6yOwNpN$11.tO5yOvMpP4.OvW5yOxJqE$12.rDwN2yOvOrGE5.rQxP4yOwSpP$14.OpAL7.qOwL4yOsK$24.rDwDyLyOvWqT!' },
    'nivium': { name: 'Nivium incarceratus', code: '3R6i', R: 18, betas: [1 / 2, 0, 1], core: 'quad4', growth: 'poly', mu: 0.2, sigma: 0.027, dt: 0.1,
      rle: '20.O$$$$$$26.sHpNJpD$25.wXyOvOrLrDB$11.qE13.3yOvTrA$10.2yOyBxJxMwL9.4yOtDpK$9.sW6yOuP6.BrIyOyIuIyOvMqOO$9.yL4yOrAvJyOrXpS2qCpUqErVuCxMtV.tQqJqTpA$9.xJ2yOyL3.tLyOvTuXuFtOtDuCuFvR2.rVrIqEE$10.uF2yO4.xEyOwIvJuKuFvW2yO2.rGuKsC$L10.OyOsR3.sM7yOxJtQrQtQyOrQ$12.uKyLqE2.rI7yOwFvEwI2yOrD$12.pN2yOxHsFrVwFyOwFpK3.wAxCyI2yOwLqH$12.TtV2yOsP9.uI4yOwXrQ$12.TtJyOxWT10.xW3yOxRwQsWqJpF2VB$11.OpSxJyOtQ11.pXyL5yOwXuSvW2yO$10.GqCqR2yO13.V4yO2yD4yO$9.GqMrGrS2yO9.G4.rQ3yOsRvMyB2yOvB$5.qWrX3sWsUuKwQyOwA9.GpD4.2yOwST.vR2yOqE$4.uN2yOwLuUvJxH2yOyI4.Q5.pStQJ3.yOyLyOvWvMyDuUrA$4.yIyOyLsFqOrLtTxCxWuI4.qRpP4.uKxWsC3.yO2wNyGvO$3.vO3yOsH2.tOwNxEsP4.sWtB3.sU2yOrN3.yOxRtVuA$3.vO4yOsHwSxU2yOwS4.uCuPtOtDwL2yOxJ3.OyOyIqJJ$3.qE10yOuI3.qJtGuSwNxWyLxRpF3.2yOvO$4.uAxWvB4.rLxH2yOtT5.sWwLvRE3.G2yOtD$13.sH2yOqW11.2yOyB$14.rQ2yO10.xW2yOxMpN$15.xRyOxC8.xM5yOE$15.tB7yOxRvT7yOrV$15.tQ13yOuFxHyIvGQ$14.G5yOqRTpKL.qHwA2yOrD.xUyOrXL$14.rQyOxMyByOxR8.xRyOuU2yOyDqRL$.pSpF10.pPvWyDsPvRyOtG9.wD3yOxWpU$.qJuSJ8.pPxJyOyGuPxUyO10.BqTtLsCqR$2.wLuP7.pFvB3yOyI2yO10.BqJrNrDT$2.qMvRqJ7.sPwDvWwF2yOsP10.BpSqET$3.pSrL8.rSuNuPuSvEsP$15.pPrGsHrN$$$$$26.2B$23.EpPqOqCL$23.EpAT!' },
    'trilapillium': { name: 'Trilapillium inversus', code: '2L3i', R: 18, betas: [1, 1 / 3], core: 'quad4', growth: 'poly', mu: 0.25, sigma: 0.03, dt: 0.1,
      rle: '9.2A$6.MqDrKsArWrFpTG$5.pTsHtRuFuIuGtUsVrCT$4.pXsSuCuKuQxF3yOxEsGW$3.UsAtUwG8yOyNpB$3.qCtH7yOyJyI3yOpB$2.CrG7yOyFxBsMuDyB2yOX$2.P2yOyMyEyH2yOyMxG3.UyHyOvWX$2.rS2yOxCrVuHxSyJxXqV4.xH2yOrMX$2.2yOxM3.PsMqNpI4.xJ2yOuJrEF$.G2yOqR10.tMyE2yOwKsTpN$.X2yOsC8.qJwHyK3yOwRtQqS$.pS2yOxNqB7.uA5yOvMuCrF$.qL2yOyJxTxOxMvP4.tQyM4yOuPuEqX$BrCxK6yOuJ4.xR4yOuOtWqD$.rGtP6yOwB4.vByH2yOxTuLsVT$.qQtPvD4yOyMxD4.vVyG2yOuJtPqV$.pEtAuLvD3yOyMxRpC3.xK2yOxLtHrNM$2.qWtXuKuI3yOyLxVwQwPyE2yOyHsFqSO$2.GrLtUuEuA8yOwAqGpDB$3.FqTsRtItGvTyA3yOuESD$5.pAqLqWqOpTVF!' },
  };
  const PATTERN_R_MIN = 12;   // below this, Heptadentium dies (numpy)
  // A pattern by its name: the key, or the full name or its first word, case and spaces ignored
  const norm = x => String(x).toLowerCase().replace(/\s+/g, '');
  function patternKey(name) {
    const n = norm(name);
    return Object.keys(PATTERNS).find(k => k === n || norm(PATTERNS[k].name) === n || norm(PATTERNS[k].name.split(' ')[0]) === n) ?? null;
  }
  const GONE = 0.4;      // a dot is eaten when its window holds less than this of its stamped mass
  const POINTS = 10;     // per dot the player's ship eats
  // Colour: a pale blue at first, its hue turning slowly towards cyan and on round the wheel.
  const HUE0 = 220, HUE_DEG_PER_S = -2, SAT = 0.6, LIGHT = 0.82;   // a full turn in 3 minutes
  // The field drifts as the star field (to the left by default; [stars] sets its angle and speed),
  // a bit slower than its near layer at speed 1 (draw.glsl starLayer 3.5 cells per second), with no
  // sway, by game time (core driftSec: the simulation's rate while playing, so the dots keep pace
  // with the ship at any speed and on any device; while paused, the speed-1 rate in real time),
  // while the stars go by real time. Sites are kept in the field's own frame; a site at (col, row) is at board
  // (col - offset[0], row - offset[1]).
  const DRIFT = 6;     // cells per game second (30 steps: 1 cell per 5 steps)
  let offset = [0, 0];

  let kind = null, dotR = 9, form = null, rad = 9;   // kind: 'dots' (round dots) or a PATTERNS key
  let on = false, target = 0, sites = [], score = 0, t0 = 0;
  let kept = [];   // [keep_beacon] beacons eaten by the ship, waiting to be put in afresh
  // The scrolling world (the game's [clear_world:scroll]; null: none): { toBoard(wx, wy, margin):
  // the board cell [col, row] of a world point, or null if it is farther than half the board less
  // margin from the ship (on either axis); toWorld(col, row); extended(reach): a world point out of
  // the world view; ok(col, row, reach): a board spot away from the seam }. A site then also has
  // its world point (wx, wy), and is parked (wiped off the board, kept by its world point) while
  // too far from the ship; it comes back when the ship nears it.
  let world = null;
  let beaconsDrift = true;   // the beacons drift with the field as the round dots do (the game: not in challenge mode)
  // cells inside the board's edge less the dot's window (about 21 cells for Heptadentium) at which a
  // site comes back / goes: parked before the seam band (6 cells, cleared each frame) reaches it, as it
  // would wear it away ("eaten"); back as near the band as that allows, mostly under the game's black
  // stripe (2 R: a beacon, 42 cells across, can't hide wholly between it and the band)
  const PARK_IN = 10, PARK_OUT = 6;
  function setWorld(w) { world = w; }
  const hooks = { score() {}, color(hue, rgb) {} };   // color: each frame while on, the dots' hue (degrees) and colour

  // The soliton resampled to `scale` of its size, bilinearly (Carlman's resizeSoliton).
  function resize(src, scale) {
    const h = src.length, w = src[0].length;
    const nh = Math.max(1, Math.round(h * scale)), nw = Math.max(1, Math.round(w * scale));
    const flat = new Float32Array(nh * nw);
    let mass = 0;
    for (let y = 0; y < nh; y++) {
      const fy = (y + 0.5) * h / nh - 0.5;
      const y0 = Math.max(0, Math.min(h - 1, Math.floor(fy))), y1 = Math.min(h - 1, y0 + 1), ty = Math.min(1, Math.max(0, fy - y0));
      for (let x = 0; x < nw; x++) {
        const fx = (x + 0.5) * w / nw - 0.5;
        const x0 = Math.max(0, Math.min(w - 1, Math.floor(fx))), x1 = Math.min(w - 1, x0 + 1), tx = Math.min(1, Math.max(0, fx - x0));
        const a = src[y0][x0] + (src[y0][x1] - src[y0][x0]) * tx, b = src[y1][x0] + (src[y1][x1] - src[y1][x0]) * tx;
        flat[y * nw + x] = a + (b - a) * ty;
        mass += flat[y * nw + x];
      }
    }
    return { flat, h: nh, w: nw, mass };
  }

  // The channel's rule and stamp: k = 'dots' (the round dots) or a PATTERNS key. A change of rule
  // empties the field (what is there would not live under the other rule).
  function setup(k) {
    if (kind === k) return;
    const changed = kind !== null;
    kind = k;
    if (k === 'dots') {
      dotR = Math.max(3, Math.round(C.R / 2));   // half the creatures' radius, as Carlman's dots
      form = resize(SOLITON, dotR / SOLITON_R);
      SimGL.setDotRule({ R: dotR, ...RULE });
    } else {
      const P = PATTERNS[k];
      if (!P.cells) P.cells = Rule.rle2cells(P.rle);
      dotR = Math.max(PATTERN_R_MIN, Math.round(P.R * C.R / 18));   // its own R at the game's R 18
      form = resize(P.cells, dotR / P.R);
      SimGL.setDotRule({ R: dotR, mu: P.mu, sigma: P.sigma, dt: P.dt, betas: P.betas, core: P.core, growth: P.growth });
    }
    rad = Math.max(dotR, Math.ceil(Math.max(form.h, form.w) / 2) + 1);   // the window holds the whole stamp
    if (changed && on) { sites = []; kept = []; SimGL.setDots(true); }
  }

  // A spot for a new dot: inside the field's frame (the stamp does not wrap), at least 3 dot radii
  // from the other dots and clear of every creature. how: random (the default), far (of 30 such
  // spots, the one farthest from the player's ship), center (the board's centre, as it is now:
  // the field drifts). Null if none was found this time.
  function findSpot(how = 'random') {
    const m = rad + 1, W = C.BW, H = C.BH, comps = C.found ? C.found.comps : [];
    const wrapF = (v, n) => ((v % n) + n) % n;
    if (how === 'center') {
      const c = world ? world.center() : [W / 2, H / 2];   // (the scrolling world: the view's middle)
      const col = Math.round(wrapF(c[0] + offset[0], W)), row = Math.round(wrapF(c[1] + offset[1], H));
      return world ? { col, row } : { col: Math.min(W - m, Math.max(m, col)), row: Math.min(H - m, Math.max(m, row)) };   // (stamps wrap now; the torus's frame kept as it was)
    }
    const ship = how === 'far' ? C.ship() : null;
    let best = null, bestD = -1, found = 0;
    for (let k = 0; k < 400 && found < (ship ? 30 : 1); k++) {
      const col = m + Math.floor(Math.random() * (W - 2 * m)), row = m + Math.floor(Math.random() * (H - 2 * m));
      if (world && !world.ok(col - offset[0], row - offset[1], rad + PARK_IN)) continue;   // away from the scrolling world's seam
      if (sites.some(s => Math.hypot(C.wrapD(s.col - col, W), C.wrapD(s.row - row, H)) < 3 * dotR)) continue;
      if (comps.some(c => Math.hypot(C.wrapD(c.col - (col - offset[0]), W), C.wrapD(c.row - (row - offset[1]), H)) < (c.radius || C.R) + C.R + rad)) continue;
      found++;
      const d = ship ? Math.hypot(C.wrapD(ship.col - (col - offset[0]), W), C.wrapD(ship.row - (row - offset[1]), H)) : 0;
      if (d > bestD) { best = { col, row }; bestD = d; }
    }
    return best;
  }
  // (points: what the ship eating it scores; default POINTS)
  // how extended: out of the world view (the scrolling world: a parked site, put on the board as
  // the ship nears it; else as far)
  function spawn(how, points = POINTS, dist) {
    if (how === 'extended') {
      if (world && kind !== 'dots') { const [wx, wy] = world.extended(rad, dist); sites.push({ wx, wy, parked: true, mass: form.mass, points }); return true; }
      how = 'far';
    }
    let p = null;
    const at = how && how.place ? how.place(form.w, form.h) : Array.isArray(how) ? how : null;   // a placer (the scrolling world's edge), or a board spot
    if (at) p = { col: Math.round((((at[0] + offset[0]) % C.BW) + C.BW) % C.BW), row: Math.round((((at[1] + offset[1]) % C.BH) + C.BH) % C.BH) };
    if (!p) p = findSpot(how);
    if (!p) return false;
    stamp(p.col, p.row);
    const s = { col: p.col, row: p.row, mass: form.mass, points };
    if (world && kind !== 'dots') [s.wx, s.wy] = world.toWorld(p.col - offset[0], p.row - offset[1]);   // (a pattern, a beacon: kept by its world point; round dots are not: the seam takes them, others are put in)
    sites.push(s);
    return true;
  }
  // Write an h x w block (flat, or zeros if null) into the field at (row0, col0), wrapping round the
  // torus (putDots takes only blocks inside the board)
  function putWrapped(flat, h, w, row0, col0) {
    const W = C.BW, H = C.BH;
    row0 = ((row0 % H) + H) % H; col0 = ((col0 % W) + W) % W;
    for (const [y0, y1] of [[0, Math.min(h, H - row0)], [Math.min(h, H - row0), h]])
      for (const [x0, x1] of [[0, Math.min(w, W - col0)], [Math.min(w, W - col0), w]]) {
        if (y1 <= y0 || x1 <= x0) continue;
        let part = null;
        if (flat) { part = new Float32Array((y1 - y0) * (x1 - x0)); for (let y = y0; y < y1; y++) for (let x = x0; x < x1; x++) part[(y - y0) * (x1 - x0) + x - x0] = flat[y * w + x]; }
        SimGL.putDots(part, y1 - y0, x1 - x0, (row0 + y0) % H, (col0 + x0) % W);
      }
  }
  const stamp = (col, row) => putWrapped(form.flat, form.h, form.w, row - Math.floor(form.h / 2), col - Math.floor(form.w / 2));
  // Wipe what is left of a dot (an eaten one's sliver, which might otherwise regrow).
  function wipe(s) {
    const e = 2 * rad + 1;
    putWrapped(null, e, e, s.row - rad, s.col - rad);
  }
  // The scrolling world, each frame: a site too far from the ship is parked, a parked one near
  // enough put back on the board (its world point kept as it drifts with the field)
  function scrollSites() {
    if (!world) return;
    const W = C.BW, H = C.BH;
    for (const s of sites) {
      if (s.wx === undefined) {   // a round dot touching an edge: taken away, its replacement queued for the opposite edge
        const q = world.edgeOut(s.col - offset[0], s.row - offset[1], rad);
        if (q) { s.gone = true; queue.push(q); }
        continue;
      }
      if (s.parked) {
        const b = world.toBoard(s.wx, s.wy, rad + PARK_IN);
        if (!b) continue;
        s.col = Math.round((((b[0] + offset[0]) % W) + W) % W); s.row = Math.round((((b[1] + offset[1]) % H) + H) % H);
        stamp(s.col, s.row); s.parked = false;
      } else {
        [s.wx, s.wy] = world.toWorld(s.col - offset[0], s.row - offset[1]);
        if (!world.toBoard(s.wx, s.wy, rad + PARK_OUT)) { wipe(s); s.parked = true; }
      }
    }
  }

  // [dot:<n>]: keep n dots (the score starts from 0 if the dots were off); [dot:off]: none, no score.
  // how: where the first n are put (findSpot's; then the replacements too, but in a scrolling world).
  // The scrolling world: a dot touching an edge is taken away (no score) and its replacement put on
  // the opposite edge (world.edgeOut, queue), other replacements (an eaten dot's) on the edge the ship
  // heads to (world.edgeSpot)
  let refill = 'random', puts = 0, queue = [];   // puts: dots put in since [dot:<n>] (the first n by how)
  function set(n, how = 'random') {
    if (n === 'off') {
      on = false; target = 0; sites = []; kept = [];
      SimGL.setDots(false);
      hooks.score();
      return;
    }
    turnOn();
    setup('dots');
    target = Math.min(n, SimGL.MAX_SITES);
    refill = how; puts = sites.length; queue = [];
    while (sites.length > target) wipe(sites.pop());
    hooks.score();
  }

  function turnOn() {
    if (on) return;
    on = true; score = 0; sites = []; kept = []; t0 = performance.now(); SimGL.setDots(true);
  }
  // [dot:<pattern>[:<n>][:<how>]]: put one <pattern> (a PATTERNS name, e.g. Heptadentium) in the dots
  // channel at a spot found by `how` (random, far, center: findSpot), switching the channel to its
  // rule (the round dots go). Eaten by the ship, it scores `points` (default as a dot); it is not
  // replaced. False if the name is unknown or no spot was found.
  // dist: (extended, a scrolling world) its distance from the ship, in R
  // keep ([keep_beacon]: the PATTERNS keys it may be): once the ship eats it, another is put in
  // afresh, by the same how (a new spot), again and again: a pattern picked at random from keep, not
  // the one just eaten (if keep has another), while no other beacon is on the board (one rule at a
  // time: it would empty the channel); else the same pattern
  function put(name, how = 'random', points = POINTS, dist, keep = null) {
    const k = patternKey(name);
    if (!k) return false;
    turnOn();
    setup(k);
    target = 0;   // nothing topped up
    const n = sites.length, ok = spawn(how, points, dist);
    if (ok && keep) sites[n].keep = { how, dist, list: keep };
    return ok;
  }

  // Each frame: score and replace the eaten dots, top up to the target, turn the hue.
  function frame() {
    if (!on) return;
    const dt = kind !== 'dots' && !beaconsDrift ? 0 : C.driftSec, d = C.starDrift();   // (beacons: stationary unless beaconsDrift)
    const cells = d.cells ?? DRIFT * d.speed;   // cells per game second ([stars:follow]: the ship's speed times its factor)
    offset = [(offset[0] - d.dir[0] * cells * dt) % C.BW, (offset[1] - d.dir[1] * cells * dt) % C.BH];
    SimGL.setDotOffset(...offset);
    scrollSites();
    limitMass();
    const live = sites.filter(s => !s.parked);
    const left = SimGL.readDotSites(live.map(s => [s.col, s.row]), rad);
    const found = C.found;
    // the ship's parts (a twin or trio ship: all of them)
    const shipParts = C.focusIdx >= 0 ? (C.objects().find(o => o.parts.includes(C.focusIdx))?.parts ?? [C.focusIdx]) : [];
    const lost = [];   // beacons gone, but not to the ship: put back nearby
    sites = sites.filter(s => {
      if (s.parked) return true;
      if (s.gone) { wipe(s); return false; }   // (not eaten: no score)
      const i = live.indexOf(s);
      if (left[i] >= GONE * s.mass) {
        if (kind !== 'dots' && left.com) {   // a beacon: its site follows the pattern's centre of mass (a moving one stays in its window)
          const [dx, dy] = left.com[i];   // (from the window's centre, the site rounded)
          s.col = (((Math.floor(s.col + 0.5) + dx) % C.BW) + C.BW) % C.BW; s.row = (((Math.floor(s.row + 0.5) + dy) % C.BH) + C.BH) % C.BH;
        }
        return true;
      }
      // eaten: by the player's ship if it is the creature there (or the nearest one close by)
      const at = found ? Creatures.pick(found, s.col - Math.round(offset[0]), s.row - Math.round(offset[1]), C.R + rad, SimGL.BW, SimGL.BH, SimGL.BLOCK, C.BW, C.BH) : -1;
      const byShip = at >= 0 && shipParts.includes(at);
      if (byShip) { score += s.points ?? POINTS; hooks.score(); if (s.keep) kept.push(s); }
      wipe(s);
      if (!byShip && kind !== 'dots') lost.push(s);   // a beacon eaten by another creature (or worn away): not counted
      return false;
    });
    for (const s of lost) respawnNear(s);
    // [keep_beacon]: one eaten by the ship comes back at a new spot (no spot found this frame: tried again the next)
    for (let k = kept.length; k-- > 0;) {
      const s = kept[k], n = sites.length;
      if (kind === 'dots') continue;
      const pool = n === 0 ? (s.keep.list || []).filter(x => x !== kind && PATTERNS[x]) : [];   // (another pattern only with no other beacon about)
      if (pool.length) { const waiting = kept; setup(pool[Math.floor(Math.random() * pool.length)]); kept = waiting; }   // (setup empties kept: the rule changed)
      if (!spawn(s.keep.how, s.points, s.keep.dist)) continue;
      sites[n].keep = s.keep; kept.splice(k, 1);
    }
    while (sites.length < target) {   // the scrolling world: the opposite edge for one taken away at an edge; once n are in, the edge ahead
      const at = world && puts >= target ? (queue.length ? world.queuedSpot(queue.shift(), rad) : world.edgeSpot(rad)) : null;
      if (!spawn(at || refill)) break;
      puts++;
    }
    const hue = (((HUE0 + HUE_DEG_PER_S * (performance.now() - t0) / 1000) % 360) + 360) % 360;
    const rgb = hsl(hue, SAT, LIGHT);
    SimGL.setDotColor(rgb);
    hooks.color(hue, rgb);
  }

  // A beacon lost (eaten by another creature, or worn away) put back near where it was: its centre
  // 2 to 4 R from there, at a random angle, clear of every creature (else as near as could be found)
  function respawnNear(s) {
    const comps = C.found ? C.found.comps : [], W = C.BW, H = C.BH;
    let best = null, bestGap = -Infinity;
    for (let k = 0; k < 60; k++) {
      const a = Math.random() * 2 * Math.PI, d = (2 + 2 * Math.random()) * C.R;
      const col = s.col + d * Math.cos(a), row = s.row + d * Math.sin(a);   // (the field's frame)
      if (world && !world.ok(col - offset[0], row - offset[1], rad + PARK_IN)) continue;   // (the scrolling world: away from its edges)
      const gap = Math.min(Infinity, ...comps.map(c => Math.hypot(C.wrapD(c.col - (col - offset[0]), W), C.wrapD(c.row - (row - offset[1]), H)) - (c.radius || C.R) - rad));
      if (gap > bestGap) { best = [col, row]; bestGap = gap; }
      if (gap > C.R) break;   // clear: take it
    }
    const [col, row] = best || [s.col, s.row];
    const t = { col: Math.round((((col % W) + W) % W)), row: Math.round((((row % H) + H) % H)), mass: form.mass, points: s.points };
    if (world) [t.wx, t.wy] = world.toWorld(t.col - offset[0], t.row - offset[1]);
    stamp(t.col, t.row);
    sites.push(t);
  }
  // Beacons left (pattern sites, on the board or not): [need:clear_beacon] waits for none
  const beaconCount = () => kind && kind !== 'dots' ? sites.length + kept.length : 0;   // (a [keep_beacon] one waiting to come back too)

  // The field's mass limit (as the bullets' [bullet_limit], with no tag of its own): round dots, n + 2
  // dots' worth ([dot:<n>]); beacons, 3 beacons' worth (each on the board). Over it (the field
  // grew, "exploded"), it is emptied and put back as it should be: round dots afresh (frame refills
  // them; no score), beacons stamped afresh at their sites (none lost, nor counted as eaten).
  function limitMass() {
    if (!form) return;
    const live = sites.filter(s => !s.parked);
    const limit = (kind === 'dots' ? target + 2 : 3 * Math.max(1, live.length)) * form.mass;
    const blk = SimGL.readDotBlocks();
    let total = 0;
    for (let i = 0; i < blk.length; i++) total += blk[i];
    if (total <= limit) return;
    console.info(`Cognoship: the dots channel grew to ${(total / form.mass).toFixed(1)} ${kind === 'dots' ? 'dots' : 'beacons'}' worth: put back as it was`);
    SimGL.setDots(true);   // (emptied)
    if (kind === 'dots') sites = sites.filter(s => s.parked);
    else for (const s of live) stamp(s.col, s.row);
  }

  // HSL (h in degrees, s and l in 0..1) to RGB in 0..1.
  function hsl(h, s, l) {
    const a = s * Math.min(l, 1 - l);
    const f = n => { const k = (n + h / 30) % 12; return l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1)); };
    return [f(0), f(8), f(4)];
  }

  // Empty the field and put the dots in afresh (the score kept): the game's Reset.
  function clear() { if (!on) return; sites = []; kept = []; SimGL.setDots(true); }
  // The parked sites' world points (the scrolling world: the game points to them)
  const parked = () => sites.filter(s => s.parked).map(s => [s.wx, s.wy]);
  // The beacons on the board (pattern sites, not parked): [col, row, radius] (board cells), for the
  // game to put creatures clear of them
  const beaconCells = () => kind && kind !== 'dots' ? sites.filter(s => !s.parked).map(s => [s.col - offset[0], s.row - offset[1], rad]) : [];
  // A beacon's radius (its window), cells
  const beaconRadius = () => rad;
  // Every pattern site's world point (the scrolling world's beacons, on the board or not)
  const beacons = () => sites.filter(s => s.wx !== undefined).map(s => [s.wx, s.wy]);

  // As if the player's ship ate a dot (?dev=1's P key): the score goes up by POINTS
  function eat() { score += POINTS; hooks.score(); }
  return { set, put, patternKey, PATTERN_KEYS: Object.keys(PATTERNS), get kind() { return kind; }, clear, frame, eat, hooks, hsl, setWorld, parked, beacons, beaconCells, beaconRadius, beaconCount,
           set beaconsDrift(v) { beaconsDrift = !!v; }, get on() { return on; }, get score() { return score; }, get count() { return sites.length; } };
})();
