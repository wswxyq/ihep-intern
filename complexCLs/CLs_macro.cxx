#include "StandardHypoTestInvDemo.C"

void CLs_macro() {

/// type = 0 Freq calculator
/// type = 1 Hybrid calculator
/// type = 2 Asymptotic calculator
/// type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)
///
/// testStatType = 0 LEP
///              = 1 Tevatron
///              = 2 Profile Likelihood two sided
///              = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)
///              = 4 Profile Likelihood signed ( pll = -pll if mu < mu_hat)
///              = 5 Max Likelihood Estimate as test statistic
///              = 6 Number of observed event as test statistic
/// ~~~

    optHTInv.noSystematics=true;
    StandardHypoTestInvDemo("ws2file.root","w","S+B_modelNM","B_modelNM","ds",
                            2, 3, true,
                            10, 0, 100, 2000, true);
}