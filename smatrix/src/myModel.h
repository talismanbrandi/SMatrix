/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#pragma warning disable 2196
#ifndef MYMODEL_H
#define	MYMODEL_H

#include <HEPfit.h>

/**
 * @class myModel
 * @brief My own Model.
 */
class myModel: public StandardModel {
public:

    static const int NmyModelvars = 22; /* Define number of mandatory parameters in the model. */
    static const std::string myModelvars[NmyModelvars]; /* Vector of model variable names. */
    
    /**
     * @brief myModel constructor
     */
    myModel();
    
    /**
     * @brief myModel destructor
     */
    ~myModel();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual bool setFlag(const std::string name, const bool value);
    
    
    /**
     * 
     * @return the coupling Api0
     */
    double getApi0() const
    {
        return Api0;
    }

    /**
     *
     * @return the coupling Api2
     */
    double getApi2() const
    {
        return Api2;
    }

    /**
     * 
     * @return the coupling AK0
     */
    double getAK0() const
    {
        return AK0;
    }

    /**
     *
     * @return the coupling A4pi0
     */
    double getA4pi0() const
    {
        return A4pi0;
    }
    
    /**
     *
     * @return the coupling Aetet_r
     */
    double getAetet_r() const
    {
        return Aetet_r;
    }
    
    /**
     *
     * @return the coupling Aetetp_r
     */
    double getAetetp_r() const
    {
        return Aetetp_r;
    }

    /**
     *
     * @return the coupling Aetet_i
     */
    double getAetet_i() const
    {
        return Aetet_i;
    }
    
    /**
     *
     * @return the coupling Aetetp_i
     */
    double getAetetp_i() const
    {
        return Aetetp_i;
    }

    /**
     *
     * @return the coupling delta1
     */
    double getdelta1() const
    {
        return delta1;
    }

    /**
     *
     * @return the coupling delta2
     */
    double getdelta2() const
    {
        return delta2;
    }

    /**
     *
     * @return the coupling n1r
     */
    double getn1r() const
    {
        return n1r;
    }

    /**
     *
     * @return the coupling n2r
     */
    double getn2r() const
    {
        return n2r;
    }

    /**
     *
     * @return the coupling n3r
     */
    double getn3r() const
    {
        return n3r;
    }

    /**
     *
     * @return the coupling n4r
     */
    double getn4r() const
    {
        return n4r;
    }

    /**
     *
     * @return the coupling n5r
     */
    double getn5r() const
    {
        return n5r;
    }


    /**
     *
     * @return the coupling n1i
     */
    double getn1i() const
    {
        return n1i;
    }

    /**
     *
     * @return the coupling n2i
     */
    double getn2i() const
    {
        return n2i;
    }

    /**
     *
     * @return the coupling n3i
     */
    double getn3i() const
    {
        return n3i;
    }

    /**
     *
     * @return the coupling n4i
     */
    double getn4i() const
    {
        return n4i;
    }

    /**
     *
     * @return the coupling n5i
     */
    double getn5i() const
    {
        return n5i;
    }

    /**
     *
     * @return the coupling nd1
     */
    double getnd1() const
    {
        return nd1;
    }

    /**
     *
     * @return the coupling nd2
     */
    double getnd2() const
    {
        return nd2;
    }
    
    /**
     *
     * @return the coupling cA
     */
    bool get_condition_flag() const
    {
        return condition;
    }


protected:
    
    virtual void setParameter(const std::string, const double&);

private:
    
    double Api0, Api2, AK0, A4pi0, Aetet_r, Aetetp_r, Aetet_i, Aetetp_i, delta1, delta2; /* Model Parameters */
    double n1r, n2r, n3r, n4r, n5r, n1i, n2i, n3i, n4i, n5i, nd1, nd2; /* Model Parameters */
    bool condition;
    
};

#endif	/* MYMODEL_H */

