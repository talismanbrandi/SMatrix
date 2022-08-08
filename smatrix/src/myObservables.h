/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#pragma warning disable 2196
#ifndef MYOBSERVABLES_H
#define	MYOBSERVABLES_H

#include <HEPfit.h>
#include "myModel.h"


/**
 * @class myObservables
 * @brief A class for the gg -> 4l.
 */
class myObservables : public ThObservable {
public:
    myObservables(const StandardModel& SM_i);
    virtual ~myObservables();
    void updateParameters();
    
protected:
    
    /* Define the parameters here. */
    double Api0, Api2, AK0, A4pi0, Aetet_r, Aetetp_r, Aetet_i, Aetetp_i, delta1, delta2; /* Model Parameters */
    double n1r, n2r, n3r, n4r, n5r, n1i, n2i, n3i, n4i, n5i, nd1, nd2; /* Model Parameters */
    gslpp::complex eIdelta_1, eIdelta_2;

    double GF;
    gslpp::complex ii;
    gslpp::matrix<gslpp::complex> VCKM;

    // CKM Parameters
    gslpp::complex V_ud, V_us;
    gslpp::complex V_cd, V_cs;
    gslpp::complex V_ub, V_cb;
    gslpp::complex ckm_scs, ckm_scs_bar;
    gslpp::complex ckm_CPV, ckm_CPV_bar;

    double MD0, MDP;
    double MPIP, MPIM, MKP, MKM;
    double MPI0, MKS, MKL;
    double MET, METP;
    double GD0, GDP;

    double sqrt_2, sqrt_3;

    double PS(const double M, const double m1, const double m2) const;
    double BR_bar(double amp, double amp_bar);
    double aCP(double amp, double amp_bar);

    gslpp::complex Api_0, AK_0, Arho_0, Aetet, Aetetp;
    gslpp::complex Api_0_bar, AK_0_bar, Arho_0_bar, Aetet_bar, Aetetp_bar;
    gslpp::complex Api_2, AK_11, Api_2_bar, AK_11_bar;
    gslpp::complex AK13, AK13_bar;

private:
    const myModel * my_model;
    gslpp::matrix<gslpp::complex> O_Matrix;
    gslpp::vector<gslpp::complex> phiR;
    gslpp::vector<gslpp::complex> phiI;
    gslpp::matrix<gslpp::complex> FSIR;
    gslpp::matrix<gslpp::complex> FSII;
    gslpp::complex r_CKM;

    gslpp::complex Api0_r, Api0_i, AK11, AK0_r, AK0_i, A4pi0_r, A4pi0_i;
    gslpp::vector<gslpp::complex> A0T_r, A0T_i;

};

class BRpppm : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRpppm(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRpppm
     */
    double computeThValue ();

private:
    int obstype;

};

class BRp0p0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRp0p0(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRp0p0
     */
    double computeThValue ();

private:
    int obstype;

};

class BRppp0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRppp0(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRppp0
     */
    double computeThValue ();

private:
    int obstype;

};

class BRkpkm : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRkpkm(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRkpkm
     */
    double computeThValue ();

private:
    int obstype;

};

class BRkSkS : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRkSkS(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRksks
     */
    double computeThValue ();

private:
    int obstype;

};

class BRkpkS : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRkpkS(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRkpkS
     */
    double computeThValue ();

private:
    int obstype;

};

class BRetet : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRetet(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRetet
     */
    double computeThValue ();

private:
    int obstype;

};

class BRetetp : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRetetp(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRetetp
     */
    double computeThValue ();

private:
    int obstype;

};


class DACP : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    DACP(const StandardModel& SM_i);

    /**
     * @return DACP
     */
    double computeThValue ();

private:

    BRpppm myBRpppm;
    BRkpkm myBRkpkm;

};

class chisq_acp : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    chisq_acp(const StandardModel& SM_i);

    /**
     * @return chisq_acp
     */
    double computeThValue ();

private:

    BRpppm myBRpppm;
    BRkpkm myBRkpkm;
    BRp0p0 myBRp0p0;
    BRkSkS myBRkSkS;
    DACP myDACP;
    double mu_acp_pppm, mu_acp_kpkm, mu_acp_p0p0, mu_acp_kSkS, mu_acp_DACP;
    double s_acp_pppm, s_acp_kpkm, s_acp_p0p0, s_acp_kSkS, s_acp_DACP;

};

class BR4pi_ppmm : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BR4pi_ppmm(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BR4pi_ppmm
     */
    double computeThValue ();

private:
    int obstype;
    double c0, c2;

};

class BR4pi_pmzz : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BR4pi_pmzz(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BR4pi_pmzz
     */
    double computeThValue ();

private:
    int obstype;
    double c0, c2;

};

class A_PI_0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    A_PI_0(const StandardModel& SM_i);

    /**
     * @return A_PI_0
     */
    double computeThValue ();

private:
};

class PHI_PI_0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    PHI_PI_0(const StandardModel& SM_i);

    /**
     * @return PHI_PI_0
     */
    double computeThValue ();

private:
};

class A_RHO_0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    A_RHO_0(const StandardModel& SM_i);

    /**
     * @return A_RHO_0
     */
    double computeThValue ();

private:
};

class PHI_RHO_0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    PHI_RHO_0(const StandardModel& SM_i);

    /**
     * @return PHI_RHO_0
     */
    double computeThValue ();

private:
};

class A_K_13 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    A_K_13(const StandardModel& SM_i);

    /**
     * @return A_K_13
     */
    double computeThValue ();

private:
};

class A_K_11 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    A_K_11(const StandardModel& SM_i);

    /**
     * @return A_K_11
     */
    double computeThValue ();

private:
};

class A_K_0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    A_K_0(const StandardModel& SM_i);

    /**
     * @return A_K_0
     */
    double computeThValue ();

private:
};

class PHI_K_0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    PHI_K_0(const StandardModel& SM_i);

    /**
     * @return PHI_K_0
     */
    double computeThValue ();

private:
};

class A_ETET : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    A_ETET(const StandardModel& SM_i);

    /**
     * @return A_ETET
     */
    double computeThValue ();

private:
};

class A_ETETP : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    A_ETETP(const StandardModel& SM_i);

    /**
     * @return A_ETETP
     */
    double computeThValue ();

private:
};

class PHI_ETET : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    PHI_ETET(const StandardModel& SM_i);

    /**
     * @return PHI_ETET
     */
    double computeThValue ();

private:
};

class PHI_ETETP : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    PHI_ETETP(const StandardModel& SM_i);

    /**
     * @return PHI_ETETP
     */
    double computeThValue ();

private:
};

#endif	/* MYOBSERVABLES_H */

