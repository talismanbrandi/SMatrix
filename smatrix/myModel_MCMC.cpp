/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 * @example myModel_MCMC.cpp
 * @example myModel.h
 * @example myModel.cpp
 * @example myObservables.h
 * @example myObservables.cpp
 * This is an example of how to add user-defined model and observables
 * and to perform a Bayesian analysis with the Markov Chain Monte Carlo.
 *
 */

#pragma warning disable 2196
#include <iostream>
#include <HEPfit.h>
#include <boost/bind.hpp>
#include "myModel.h"
#include "myObservables.h"

/* Necessary if MPI support is enabled during compilation. */
#ifdef _MPI
#include <mpi.h>
#endif

int main(int argc, char** argv) 
{

/* Necessary if MPI support is enabled during compilation. */
#ifdef _MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    /* In our MPI implementation the process with rank 0 is the master. */
    int rank = 0;
#endif
    
    try {
        
        if(argc != 3){
            /* Print usage and exit. */
            if (rank == 0) std::cout << "\nusage: " << argv[0] << " ModelConf.conf MonteCarlo.conf\n" << std::endl;
            return EXIT_SUCCESS;
        }

        /* Define the model configuration file.                                */
        /* Here it is passed as the first argument to the executable. The      */
        /* model configuration file provides the values with errors for the    */
        /* mandatory model parameters, as well as the list of observables,     */
        /* observables2D, correlated Gaussian observables.                     */
        /* See documentation for details.                                      */
        std::string ModelConf = argv[1];
        
        /* Define the Monte Carlo configuration file.                         */
        /* Here it is passed as the second argument to the executable. The    */
        /* Monte Carlo configuration file provides the parameters used in the */
        /* Monte Carlo run. See documentation for details.                    */
        std::string MCMCConf = argv[2];
        
        /* Define the ROOT output file (w/o extension, empty string will set it to MCout) */
        std::string FileOut = "";        
        
        /* Define the optional job tag. */
        std::string JobTag = "";
        
        /* Create objects of the classes ModelFactory and ThObsFactory */
        ThObsFactory ThObsF;
        ModelFactory ModelF;
        myModel my_model;

        /* register user-defined model named ModelName defined in class ModelClass using the following syntax: */
        ModelF.addModelToFactory("myModel", boost::factory<myModel*>() );
        
        /* register user-defined ThObservable named ThObsName defined in class ThObsClass using the following syntax: */
        ThObsF.addObsToFactory("BRpppm", boost::bind(boost::factory<BRpppm*>(), _1, 1) );
        ThObsF.addObsToFactory("BRp0p0", boost::bind(boost::factory<BRp0p0*>(), _1, 1) );
        ThObsF.addObsToFactory("BRppp0", boost::bind(boost::factory<BRppp0*>(), _1, 1) );
        ThObsF.addObsToFactory("BRkpkm", boost::bind(boost::factory<BRkpkm*>(), _1, 1) );
        ThObsF.addObsToFactory("BRkSkS", boost::bind(boost::factory<BRkSkS*>(), _1, 1) );
        ThObsF.addObsToFactory("BRkpkS", boost::bind(boost::factory<BRkpkS*>(), _1, 1) );
        ThObsF.addObsToFactory("BRkpkL", boost::bind(boost::factory<BRkpkS*>(), _1, 1) ); /*kpKL is given by KpKS*/
        ThObsF.addObsToFactory("BRetet", boost::bind(boost::factory<BRetet*>(), _1, 1) );
        ThObsF.addObsToFactory("BRetetp", boost::bind(boost::factory<BRetetp*>(), _1, 1) );

        ThObsF.addObsToFactory("BR4pi_ppmm", boost::bind(boost::factory<BR4pi_ppmm*>(), _1, 1) );
        ThObsF.addObsToFactory("BR4pi_pmzz", boost::bind(boost::factory<BR4pi_pmzz*>(), _1, 1) );

        ThObsF.addObsToFactory("aCPpppm", boost::bind(boost::factory<BRpppm*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPp0p0", boost::bind(boost::factory<BRp0p0*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPkpkm", boost::bind(boost::factory<BRkpkm*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPkSkS", boost::bind(boost::factory<BRkSkS*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPkpkS", boost::bind(boost::factory<BRkpkS*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPetet", boost::bind(boost::factory<BRetet*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPetetp", boost::bind(boost::factory<BRetetp*>(), _1, 2) );

        ThObsF.addObsToFactory("DACP", boost::bind(boost::factory<DACP*>(), _1) );

        ThObsF.addObsToFactory("chisq_acp", boost::bind(boost::factory<chisq_acp*>(), _1) );
        
        ThObsF.addObsToFactory("A_PI_0", boost::bind(boost::factory<A_PI_0*>(), _1) );
        ThObsF.addObsToFactory("PHI_PI_0", boost::bind(boost::factory<PHI_PI_0*>(), _1) );
        ThObsF.addObsToFactory("A_RHO_0", boost::bind(boost::factory<A_RHO_0*>(), _1) );
        ThObsF.addObsToFactory("PHI_RHO_0", boost::bind(boost::factory<PHI_RHO_0*>(), _1) );
        ThObsF.addObsToFactory("A_K_13", boost::bind(boost::factory<A_K_13*>(), _1) );
        ThObsF.addObsToFactory("A_K_11", boost::bind(boost::factory<A_K_11*>(), _1) );
        ThObsF.addObsToFactory("A_K_0", boost::bind(boost::factory<A_K_0*>(), _1) );
        ThObsF.addObsToFactory("PHI_K_0", boost::bind(boost::factory<PHI_K_0*>(), _1) );
        ThObsF.addObsToFactory("A_ETET", boost::bind(boost::factory<A_ETET*>(), _1) );
        ThObsF.addObsToFactory("A_ETETP", boost::bind(boost::factory<A_ETETP*>(), _1) );
        ThObsF.addObsToFactory("PHI_ETET", boost::bind(boost::factory<PHI_ETET*>(), _1) );
        ThObsF.addObsToFactory("PHI_ETETP", boost::bind(boost::factory<PHI_ETETP*>(), _1) );
        
        /* Create an object of the class MonteCarlo. */
        MonteCarlo MC(ModelF, ThObsF, ModelConf, MCMCConf, FileOut, JobTag);

        /* Do a test run if you wish to see the values of the observables      */
        /* and the correlated Gaussian observables defined in the model        */
        /* configuration file computed with the central value of the mandatory */
        /* parameters defined in the same file.                                */
        if (MCMCConf.compare("--noMC") == 0) MC.TestRun(rank);


        /* Initiate the Mote Carlo run. */
        else MC.Run(rank);

/* Necessary if MPI support is enabled during compilation. */
#ifdef _MPI
        MPI_Finalize();
#endif
        
        return EXIT_SUCCESS;
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
