/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "myObservables.h"


myObservables::myObservables(const StandardModel& SM_i)
:   ThObservable(SM_i), VCKM(3,3,0.), 
    my_model(static_cast<const myModel*> (&SM_i)),
    O_Matrix(5, 5, 0.),
    phiR(5, 0.),
    phiI(5, 0.),
    FSIR(5, 5, 0.),
    FSII(5, 5, 0.),
    A0T_r(5, 0.),
    A0T_i(5, 0.)
{
    // The Rotaion Matrix
    O_Matrix.assign(0, 0, 0.817745);
    O_Matrix.assign(0, 1, 0.430377);
    O_Matrix.assign(0, 2, -0.113613);
    O_Matrix.assign(0, 3, 0.230412);
    O_Matrix.assign(0, 4, 0.282967);
    O_Matrix.assign(1, 0, 0.0424431);
    O_Matrix.assign(1, 1, -0.524008);
    O_Matrix.assign(1, 2, 0.239839);
    O_Matrix.assign(1, 3, -0.0537274);
    O_Matrix.assign(1, 4, 0.814374);
    O_Matrix.assign(2, 0, 0.131331);
    O_Matrix.assign(2, 1, -0.491756);
    O_Matrix.assign(2, 2, 0.0347818);
    O_Matrix.assign(2, 3, 0.813265);
    O_Matrix.assign(2, 4, -0.279854);
    O_Matrix.assign(3, 0, 0.558553);
    O_Matrix.assign(3, 1, -0.482087);
    O_Matrix.assign(3, 2, 0.112211);
    O_Matrix.assign(3, 3, -0.526586);
    O_Matrix.assign(3, 4, -0.407096);
    O_Matrix.assign(4, 0, 0.0161804);
    O_Matrix.assign(4, 1, 0.256828);
    O_Matrix.assign(4, 2, 0.956958);
    O_Matrix.assign(4, 3, 0.0730079);
    O_Matrix.assign(4, 4, -0.112603);

    // The masses
    MD0 = 1.86484;
    MDP = 1.86961;
    MPIP = 0.13957018;
    MPIM = 0.13957018;
    MKP = 0.493677;
    MKM = 0.493677;
    MPI0 = 0.1349766;
    MKS = 0.497614;
    MKL = 0.497614;
    MET = 0.547862;
    METP = 0.95778;

    // The decay widths
    GD0 = 1./0.4101 * 6.58211928e-13;
    GDP = 1./1.040 * 6.58211928e-13;

    sqrt_2 = sqrt(2.);
    sqrt_3 = sqrt(3.); 
}

myObservables::~myObservables()
{}

void myObservables::updateParameters()
{
    // Calling the model parameters from the model
    Api0 = my_model->getApi0();
    Api2 = my_model->getApi2();
    AK0 = my_model->getAK0();
    A4pi0 = my_model->getA4pi0();
    Aetet_r = my_model->getAetet() * cos(my_model->getdelta_etet());
    Aetetp_r = my_model->getAetetp() * cos(my_model->getdelta_etetp());
    Aetet_i = my_model->getAetet() * sin(my_model->getdelta_etet());
    Aetetp_i = my_model->getAetetp() * sin(my_model->getdelta_etetp());
    delta1 = my_model->getdelta1();
    delta2 = my_model->getdelta2();

    // The degeneracy counts (real and imaginary)
    n1r = my_model->getn1r();
    n2r = my_model->getn2r();
    n3r = my_model->getn3r();
    n4r = my_model->getn4r();
    n5r = my_model->getn5r();

    n1i = my_model->getn1i();
    n2i = my_model->getn2i();
    n3i = my_model->getn3i();
    n4i = my_model->getn4i();
    n5i = my_model->getn5i();

    nd1 = my_model->getnd1();
    nd2 = my_model->getnd2();

    // Definition of the phases (real and imaginary)
    phiR.assign(0, gslpp::complex(1., -0.582523 + n1r * M_PI, 1));
    phiR.assign(1, gslpp::complex(1., -0.0240778 + n2r * M_PI, 1));
    phiR.assign(2, gslpp::complex(1., 0.00899414 + n3r * M_PI, 1));
    phiR.assign(3, gslpp::complex(1., -0.361503 + n4r * M_PI, 1));
    phiR.assign(4, gslpp::complex(1., -1.28638 + n5r * M_PI, 1));

    phiI.assign(0, gslpp::complex(1., -0.582523 + n1r * M_PI, 1));
    phiI.assign(1, gslpp::complex(1., -0.0240778 + n2r * M_PI, 1));
    phiI.assign(2, gslpp::complex(1., 0.00899414 + n3r * M_PI, 1));
    phiI.assign(3, gslpp::complex(1., -0.361503 + n4r * M_PI, 1));
    phiI.assign(4, gslpp::complex(1., -1.28638 + n5r * M_PI, 1));

    // The FSI matrix
    FSIR = O_Matrix.transpose() * gslpp::matrix<gslpp::complex>(phiR) * O_Matrix;
    FSII = O_Matrix.transpose() * gslpp::matrix<gslpp::complex>(phiI) * O_Matrix;

    ii = gslpp::complex(0., 1., 0);

    // The exponential of the phases
    eIdelta_1 = exp(ii * delta1);
    eIdelta_2 = exp(ii * delta2);

    // Calling the Fermi constant and the CKM Matrix from the SM
    GF = SM.getGF();
    VCKM = SM.getVCKM();

    V_ud = VCKM(0,0);
    V_us = VCKM(0,1);
    V_cd = VCKM(1,0);
    V_cs = VCKM(1,1);
    V_ub = VCKM(0,2);
    V_cb = VCKM(1,2);

    ckm_scs = (V_us * (V_cs).conjugate() - V_ud * (V_cd).conjugate());
    ckm_scs_bar = (ckm_scs).conjugate();

    ckm_CPV = (V_us * (V_cs).conjugate() + V_ud * (V_cd).conjugate());
    ckm_CPV_bar = ckm_CPV.conjugate();

    r_CKM = (V_cb.conjugate()*V_ub/ckm_scs.real()).imag();

    // Building the ununitarized amplitudes
    Api0_r = Api0;
    Api0_i = r_CKM * Api0;

    AK13 = 0.;
    AK13_bar = 0.;
    AK11 = AK0;
    AK0_r = -AK0;
    AK0_i = r_CKM * AK0;

    A4pi0_r = A4pi0;
    A4pi0_i = r_CKM * A4pi0;

    // 5 vectors for the ununitarized amplitudes
    A0T_r.assign(0, Api0_r);
    A0T_r.assign(1, AK0_r);
    A0T_r.assign(2, A4pi0_r);
    A0T_r.assign(3, Aetet_r);
    A0T_r.assign(4, Aetetp_r);

    A0T_i.assign(0, Api0_i);
    A0T_i.assign(1, AK0_i);
    A0T_i.assign(2, A4pi0_i);
    A0T_i.assign(3, Aetet_i);
    A0T_i.assign(4, Aetetp_i);

    // final state interactions
    A0T_r = FSIR * A0T_r;
    A0T_i = FSII * A0T_i;

    // final state unitarized amplitudes I = 0
    Api_0 = A0T_r(0) + ii * A0T_i(0);
    AK_0 = A0T_r(1) + ii * A0T_i(1);
    Arho_0 = A0T_r(2) + ii * A0T_i(2);
    Aetet = A0T_r(3) + ii * A0T_i(3);
    Aetetp = A0T_r(4) + ii * A0T_i(4);

    Api_0_bar = A0T_r(0) - ii * A0T_i(0);
    AK_0_bar = A0T_r(1) - ii * A0T_i(1);
    Arho_0_bar = A0T_r(2) - ii * A0T_i(2);
    Aetet_bar = A0T_r(3) - ii * A0T_i(3);
    Aetetp_bar = A0T_r(4) - ii * A0T_i(4);

    // final state unitarized amplitudes I = 2
    Api_2 = Api2 * gslpp::complex(1., delta2, true) * (1. + ii * r_CKM * gslpp::complex(1., nd2*M_PI, true));
    AK_11 = AK11 * gslpp::complex(1., delta1, true) * (1. - ii * r_CKM * gslpp::complex(1., nd1*M_PI, true));

    Api_2_bar = Api2 * gslpp::complex(1., delta2, true) * (1. - ii * r_CKM * gslpp::complex(1., nd2*M_PI, true));
    AK_11_bar = AK11 * gslpp::complex(1., delta1, true) * (1. + ii * r_CKM * gslpp::complex(1., nd1*M_PI, true));

}

double myObservables::PS(const double M, const double m1, const double m2) const
{
    double pc = sqrt( (M*M - (m1 + m2)*(m1 + m2))*(M*M - (m1 - m2)*(m1 - m2)) )/2.0/M;
    return ( pc/(8.*M_PI*M*M)*GF*GF/2.0 );
}

double myObservables::BR_bar(double amp, double amp_bar){
    return (amp + amp_bar)/2.;
}

double myObservables::aCP(double amp, double amp_bar){
    return (amp - amp_bar)/(amp + amp_bar) * 1.e2;
}

/*******************************************************************************
 * Observables: DCS two body decays                                                                *
 * ****************************************************************************/

BRpppm::BRpppm(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRpppm::computeThValue()
{
    updateParameters();

    gslpp::complex Apppm = (Api_2 + sqrt_2 * Api_0)/2./sqrt_3;
    gslpp::complex ApppmB = (Api_2_bar + sqrt_2 * Api_0_bar)/2./sqrt_3;

    if (obstype == 1) return PS(MD0,MPIP,MPIM) * BR_bar((Apppm).abs2(), (ApppmB).abs2()) / GD0;
    else if (obstype == 2) return aCP((Apppm).abs2(), (ApppmB).abs2());
    else return 0.;
}

BRp0p0::BRp0p0(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRp0p0::computeThValue()
{
    updateParameters();

    gslpp::complex Ap0p0 = (Api_2 * sqrt_2 - Api_0)/sqrt_2/sqrt_3;
    gslpp::complex Ap0p0B = (Api_2_bar * sqrt_2 - Api_0_bar)/sqrt_2/sqrt_3;

    if (obstype == 1) return PS(MD0,MPI0,MPI0) * BR_bar((Ap0p0).abs2(), (Ap0p0B).abs2())/2. / GD0;
    else if (obstype == 2) return aCP((Ap0p0).abs2(), (Ap0p0B).abs2());
    else return 0.;
}

BRppp0::BRppp0(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRppp0::computeThValue()
{
    updateParameters();

    gslpp::complex Appp0 = sqrt_3/sqrt_2 * Api_2/2.;
    gslpp::complex Appp0B = sqrt_3/sqrt_2 * Api_2_bar/2.;

    if (obstype == 1) return PS(MDP,MPIP,MPI0) * BR_bar((Appp0).abs2(), (Appp0B).abs2()) / GDP;
    else return 0.;
}

BRkpkm::BRkpkm(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRkpkm::computeThValue()
{
    updateParameters();

    gslpp::complex Akpkm = (AK13 + AK_11 + AK_0)/2.;
    gslpp::complex AkpkmB = (AK13_bar + AK_11_bar + AK_0_bar)/2.;

    if (obstype == 1) return PS(MD0,MKP,MKM) * BR_bar((Akpkm).abs2(), (AkpkmB).abs2()) / GD0;
    else if (obstype == 2) return aCP((Akpkm).abs2(), (AkpkmB).abs2());
    else return 0.;
}

BRkSkS::BRkSkS(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRkSkS::computeThValue()
{
    updateParameters();

    gslpp::complex Ak0k0bar = (-AK13 - AK_11 + AK_0)/2.;
    gslpp::complex Ak0k0barB = (-AK13_bar - AK_11_bar + AK_0_bar)/2.;

    if (obstype == 1) return PS(MD0,MKS,MKS) * BR_bar((Ak0k0bar).abs2(), (Ak0k0barB).abs2())/2. / GD0;
    else if (obstype == 2) return aCP((Ak0k0bar).abs2(), (Ak0k0barB).abs2());
    else return 0.;
}

BRkpkS::BRkpkS(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRkpkS::computeThValue()
{
    updateParameters();

    gslpp::complex Akpk0bar = -AK13/2. + AK_11;
    gslpp::complex Akpk0barB = -AK13_bar/2. + AK_11_bar;

    if (obstype == 1) return PS(MDP,MKP,MKS) * BR_bar((Akpk0bar).abs2(), (Akpk0barB).abs2())/2. / GDP;
    else if (obstype == 2) return aCP((Akpk0bar).abs2(), (Akpk0barB).abs2());
    else return 0.;
}

BRetet::BRetet(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRetet::computeThValue()
{
    updateParameters();

    // gslpp::complex Aetet = 0.;
    // gslpp::complex AetetB = ckm_scs_bar * Aetet;
    // Aetet = ckm_scs * Aetet;

    if (obstype == 1) return PS(MD0,MET,MET) * BR_bar((Aetet).abs2(), (Aetet_bar).abs2())/2. / GD0;
    else return 0.;
}

BRetetp::BRetetp(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRetetp::computeThValue()
{
    updateParameters();

    // gslpp::complex Aetetp = 0.;
    // gslpp::complex AetetpB = ckm_scs_bar * Aetetp;
    // Aetetp = ckm_scs * Aetetp;

    if (obstype == 1) return PS(MD0,MET,METP) * BR_bar((Aetetp).abs2(), (Aetetp_bar).abs2())/2. / GD0;
    else return 0.;
}

DACP::DACP(const StandardModel& SM_i)
: myObservables(SM_i), myBRpppm(SM_i, 2), myBRkpkm(SM_i, 2)
{}

double DACP::computeThValue()
{
    updateParameters();

    return (myBRkpkm.computeThValue() - myBRpppm.computeThValue());

}


/*******************************************************************************
 * Observables: DCS four body decays                                                                *
 * ****************************************************************************/

BR4pi_ppmm::BR4pi_ppmm(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   
    obstype = obstype_i; 
    c0 = 10.*sqrt(15.)/23.;
    c2 = -2.*sqrt(21.)/23.;

}

double BR4pi_ppmm::computeThValue()
{
    updateParameters();

    gslpp::complex A4pi_ppmm = /* c2*A2r +*/ c0 * Arho_0;
    gslpp::complex A4pi_ppmmB = /* c2*A2rbar +*/ c0 * Arho_0_bar;

    if (obstype == 1) return PS(MD0,MPIP,MPIM) * BR_bar((A4pi_ppmm).abs2(), (A4pi_ppmmB).abs2()) / GD0;
    else if (obstype == 2) return aCP((A4pi_ppmm).abs2(), (A4pi_ppmmB).abs2());
    else return 0.;
}


BR4pi_pmzz::BR4pi_pmzz(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   
    obstype = obstype_i; 
    c0 = sqrt(30.)/23.;
    c2 = -5.*sqrt(21./2.)/23.;

}

double BR4pi_pmzz::computeThValue()
{
    updateParameters();

    gslpp::complex A4pi_pmzz = /* c2*A2r +*/ c0 * Arho_0;
    gslpp::complex A4pi_pmzzB = /* c2*A2rbar +*/ c0 * Arho_0_bar;

    if (obstype == 1) return PS(MD0,MPIP,MPIM) * BR_bar((A4pi_pmzz).abs2(), (A4pi_pmzzB).abs2()) / GD0;
    else if (obstype == 2) return aCP((A4pi_pmzz).abs2(), (A4pi_pmzzB).abs2());
    else return 0.;
}

/*******************************************************************************
 * Observables: Amplitudes and phases                                                               *
 * ****************************************************************************/

A_PI_0::A_PI_0(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double A_PI_0::computeThValue()
 {
    updateParameters();
    return Api_0.abs();
 }

PHI_PI_0::PHI_PI_0(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double PHI_PI_0::computeThValue()
 {
    updateParameters();
    return Api_0.arg()-Api_2.arg();
 }

 A_RHO_0::A_RHO_0(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double A_RHO_0::computeThValue()
 {
    updateParameters();
    return Arho_0.abs();
 }

PHI_RHO_0::PHI_RHO_0(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double PHI_RHO_0::computeThValue()
 {
    updateParameters();
    return Arho_0.arg()/* -A2r.arg()*/;
 }

A_K_13::A_K_13(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double A_K_13::computeThValue()
 {
    updateParameters();
    return AK13.abs();
 }

A_K_11::A_K_11(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double A_K_11::computeThValue()
 {
    updateParameters();
    return AK_11.abs();
 }

A_K_0::A_K_0(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double A_K_0::computeThValue()
 {
    updateParameters();
    return AK_0.abs();
 }

PHI_K_0::PHI_K_0(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double PHI_K_0::computeThValue()
 {
    updateParameters();
    return AK_0.arg() - AK_11.arg();
 }

A_ETET::A_ETET(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double A_ETET::computeThValue()
 {
    updateParameters();
    return Aetet.abs();
 }

A_ETETP::A_ETETP(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double A_ETETP::computeThValue()
 {
    updateParameters();
    return Aetetp.abs();
 }

 PHI_ETET::PHI_ETET(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double PHI_ETET::computeThValue()
 {
    updateParameters();
    return Aetet.arg();
 }

PHI_ETETP::PHI_ETETP(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double PHI_ETETP::computeThValue()
 {
    updateParameters();
    return Aetetp.arg();
 }