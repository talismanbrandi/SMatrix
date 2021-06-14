/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "myModel.h"

/* Define mandatory model parameters here. */
const std::string myModel::myModelvars[NmyModelvars] = {"Api0", "Api2", "AK0", "A4pi0", "Aetet_r", "Aetetp_r", 
                                                        "Aetet_i", "Aetetp_i", "delta1", "delta2",
                                                        "n1r", "n1r", "n3r", "n4r", "n5r", "n1i", "n2i", "n3i", "n4i", "n5i",
                                                        "nd1", "nd2"};

myModel::myModel()
:   StandardModel()
{
    /* Define all the parameters here and port them as observables too */
    ModelParamMap.insert(std::make_pair("Api0", std::cref(Api0)));
    ModelParamMap.insert(std::make_pair("Api2", std::cref(Api2)));
    ModelParamMap.insert(std::make_pair("AK0", std::cref(AK0)));
    ModelParamMap.insert(std::make_pair("A4pi0", std::cref(A4pi0)));
    ModelParamMap.insert(std::make_pair("Aetet_r", std::cref(Aetet_r)));
    ModelParamMap.insert(std::make_pair("Aetetp_r", std::cref(Aetetp_r)));
    ModelParamMap.insert(std::make_pair("Aetet_i", std::cref(Aetet_i)));
    ModelParamMap.insert(std::make_pair("Aetetp_i", std::cref(Aetetp_i)));
    ModelParamMap.insert(std::make_pair("delta1", std::cref(delta1)));
    ModelParamMap.insert(std::make_pair("delta2", std::cref(delta2)));
    ModelParamMap.insert(std::make_pair("n1r", std::cref(n1r)));
    ModelParamMap.insert(std::make_pair("n2r", std::cref(n2r)));
    ModelParamMap.insert(std::make_pair("n3r", std::cref(n3r)));
    ModelParamMap.insert(std::make_pair("n4r", std::cref(n4r)));
    ModelParamMap.insert(std::make_pair("n5r", std::cref(n5r)));
    ModelParamMap.insert(std::make_pair("n1i", std::cref(n1i)));
    ModelParamMap.insert(std::make_pair("n2i", std::cref(n2i)));
    ModelParamMap.insert(std::make_pair("n3i", std::cref(n3i)));
    ModelParamMap.insert(std::make_pair("n4i", std::cref(n4i)));
    ModelParamMap.insert(std::make_pair("n5i", std::cref(n5i)));
    ModelParamMap.insert(std::make_pair("nd1", std::cref(nd1)));
    ModelParamMap.insert(std::make_pair("nd2", std::cref(nd2)));
}

myModel::~myModel()
{
    if (IsModelInitialized()) {
        /* Destroy whatever you want, e.g. potentially dangerous pointers. */
    }
}

/* Initialize model here */
bool myModel::InitializeModel()
{
    condition = false;
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}
    
bool myModel::Init(const std::map<std::string, double>& DPars)
{
    return(StandardModel::Init(DPars));
}

/* Do whatever is necessary before parameters are updated by the MCMC. */
bool myModel::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);
    return (true);
}

/* Model update method used be the MCMC to update the model parameters. */
bool myModel::Update(const std::map<std::string, double>& DPars)
{    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);

    return (true);
}

/* Postupdate method to update whatever is needed after the model parameters are updated */
bool myModel::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    return (true);
}

/* Model parameters and their derived quantities can be set here. */
void myModel::setParameter(const std::string name, const double& value)
{
    if(name.compare("Api0") == 0)
        Api0 = value;
    else if(name.compare("Api2") == 0)
        Api2 = value;
    else if(name.compare("AK0") == 0)
        AK0 = value;
    else if(name.compare("A4pi0") == 0)
        A4pi0 = value;
    else if(name.compare("Aetet_r") == 0)
        Aetet_r = value;
    else if(name.compare("Aetetp_r") == 0)
        Aetetp_r = value;
    else if(name.compare("Aetet_i") == 0)
        Aetet_i = value;
    else if(name.compare("Aetetp_i") == 0)
        Aetetp_i = value;
    else if(name.compare("delta1") == 0)
        delta1 = value;
    else if(name.compare("delta2") == 0)
        delta2 = value;
    else if(name.compare("n1r") == 0)
        n1r = value;
    else if(name.compare("n2r") == 0)
        n2r = value;
    else if(name.compare("n3r") == 0)
        n3r = value;
    else if(name.compare("n4r") == 0)
        n4r = value;
    else if(name.compare("n5r") == 0)
        n5r = value;
    else if(name.compare("n1i") == 0)
        n1i = value;
    else if(name.compare("n2i") == 0)
        n2i = value;
    else if(name.compare("n3i") == 0)
        n3i = value;
    else if(name.compare("n4i") == 0)
        n4i = value;
    else if(name.compare("n5i") == 0)
        n5i = value;
    else if(name.compare("nd1") == 0)
        nd1 = value;
    else if(name.compare("nd2") == 0)
        nd2 = value;    

    else
        StandardModel::setParameter(name,value);
}

bool myModel::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NmyModelvars; i++) {
        if (DPars.find(myModelvars[i]) == DPars.end()) {
            std::cout << "missing mandatory myModel parameter " << myModelvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}


/* Model Flags can be set here. */
bool myModel::setFlag(const std::string name, const bool value)
{
    bool res = false;
    
    if(name.compare("condition") == 0){
        condition = value;
        res = true;
    } else {
        res = StandardModel::setFlag(name,value);
    }

    return(res);
}
