/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTSANWITHFUNNYCURRENTBLOCK_HPP_
#define TESTSANWITHFUNNYCURRENTBLOCK_HPP_

/*
 * [[PageOutline]]
 *
 * = Sino-atrial node simulation with funny current block =
 *
 * This is the code that was used to perform the simulation in Capel ''et al.'' (2014).
 *
 * == Code Walkthrough ==
 */

/* This include is needed to utilise the cxx test framework, which we use to execute programs */
#include <cxxtest/TestSuite.h>

/* Header files from core Chaste that we will use */
#include "OutputFileHandler.hpp"
#include "Warnings.hpp"
#include "ZeroStimulus.hpp"
#include "CellProperties.hpp"
#include "dokos_model_1996Cvode.hpp"

/* All of the bulk of the simulation is performed inside this testing class */
class TestSanWithFunnyCurrentBlock : public CxxTest::TestSuite
{
private:
    /*
     * We define a helper method to:
     *  1. Run the model to a 'steady state' for this level of block (100 second run).
     *  1. Do a fine resolution run for analysis (2 second simulation)
     *  1. Calculate the cycle length and APD50 and return these to the main program below.
     */
    c_vector<double,2> RunToSteadyStateGetCycleLengthAndApd(const std::string& rExperimentName,
                                                            boost::shared_ptr<AbstractCvodeCell> pModel,
                                                            unsigned blockLevel,
                                                            unsigned experimentIndex)
    {
        std::stringstream filename;
        std::stringstream foldername;
        filename << pModel->GetSystemName() << "_block_" << blockLevel;
        foldername << rExperimentName << "_" << experimentIndex;
        unsigned voltage_index = pModel->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");

        // When solving for a long time with no output we need to increase the default number of steps CVODE can take.
        pModel->SetMaxSteps(1e6);
        pModel->Solve(0,100000,10); // Solve for t=0 to t=100 seconds, maximum time step 10ms (will be auto-refined by CVODE).

        // Get output every 0.1ms over a 2 second run for analysing carefully.
        OdeSolution solution = pModel->Solve(0,2000,0.1,0.1);
        std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);
        CellProperties voltage_properties(voltages, solution.rGetTimes(), -30);

        // Output these voltage traces to file for plotting with matlab script later
        OutputFileHandler handler(foldername.str(),false);
        out_stream p_file = handler.OpenOutputFile(filename.str());
        assert(voltages.size()==solution.rGetTimes().size());
        for (unsigned i=0; i<solution.rGetTimes().size(); i++)
        {
            *p_file << solution.rGetTimes()[i] << "\t" << voltages[i] << std::endl;
        }
        p_file->close();

        // Calculate APD50
        double apd70;
        try
        {
            apd70 = voltage_properties.GetLastActionPotentialDuration(50);
        }
        catch (Exception &e)
        {
            apd70 = -1;
        }

        // Calculate cycle length.
        std::vector<double> upstroke_times;
        double cycle_length;
        try
        {
            upstroke_times = voltage_properties.GetTimesAtMaxUpstrokeVelocity();
        }
        catch (Exception &e)
        {
            WARNING(e.GetShortMessage());
        }

        if (upstroke_times.size() > 1u)
        {
            for (unsigned i=0; i<upstroke_times.size()-1u; i++)
            {
                cycle_length = upstroke_times[i+1] - upstroke_times[i];
            }
        }
        else
        {
            cycle_length = -1;
            WARNING("Only one upstroke detected.");
        }

        c_vector<double,2> results;
        results[0] = apd70;
        results[1] = cycle_length;
        return results;
    }

    /* Define the main function that will be run using the testing framework */
public:
    void TestSanAction() throw (Exception)
    {
        // Set up a SAN model with a CVODE solver.
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        // Some models will accept a stimulus so we explicitly make it zero.
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new ZeroStimulus);
        // Create a cell model with a CVODE solver
        boost::shared_ptr<AbstractCvodeCell> p_model(new Celldokos_model_1996FromCellMLCvode(p_solver, p_stimulus));

        std::cout << "\nModel = " << p_model->GetSystemName() << "\n";

        // Get the control funny current conductances for sodium and potassium, IK and ICaL.
        const double gf_Na = p_model->GetParameter("membrane_hyperpolarisation_activated_funny_current_sodium_component_conductance");
        const double gf_K = p_model->GetParameter("membrane_hyperpolarisation_activated_funny_current_potassium_component_conductance");
        const double g_K = p_model->GetParameter("membrane_delayed_rectifier_potassium_current_conductance");
        const double g_CaL = p_model->GetParameter("membrane_L_type_calcium_current_conductance");

        // Now we'll do four experiments for fun - block of each current, and then block of them in the correct proportions.
        for (unsigned channels=0; channels<3; channels++)
        {
            std::cout << "\nNew experiment: " << channels  << "\n";
            // Reset to defaults
            p_model->SetParameter("membrane_hyperpolarisation_activated_funny_current_sodium_component_conductance", gf_Na);
            p_model->SetParameter("membrane_hyperpolarisation_activated_funny_current_potassium_component_conductance", gf_K);
            p_model->SetParameter("membrane_delayed_rectifier_potassium_current_conductance", g_K);
            p_model->SetParameter("membrane_L_type_calcium_current_conductance", g_CaL);

            // Loop over different scaling factors, and set the conductance parameters accordingly.
            for (double scaling = 1.0; scaling>=0.0; scaling -= 0.1)
            {
                if (channels==0)
                {
                    p_model->SetParameter("membrane_hyperpolarisation_activated_funny_current_sodium_component_conductance", gf_Na*scaling);
                    p_model->SetParameter("membrane_hyperpolarisation_activated_funny_current_potassium_component_conductance", gf_K*scaling);
                }
                else if (channels==1)
                {
                    p_model->SetParameter("membrane_delayed_rectifier_potassium_current_conductance", g_K*scaling);
                }
                else if (channels==2)
                {
                    p_model->SetParameter("membrane_L_type_calcium_current_conductance", g_CaL*scaling);
                }
                else
                {
                    EXCEPTION("Not implemented");
                }

                unsigned block_level = round(100*(1.0-scaling));
                c_vector<double,2> results = RunToSteadyStateGetCycleLengthAndApd("SanWithFunnyCurrentBlock", p_model, block_level, channels);

                std::cout << "Block = " << block_level << "%,\tCycle length = " << results[1] << "ms,\tAPD50 = " << results[0] << "ms\n";
            }
        }
    }

    void TestActionAt3uM() throw (Exception)
    {
        // Set up a SAN model with a CVODE solver.
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        // Some models will accept a stimulus so we explicitly make it zero.
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new ZeroStimulus);
        // Create a cell model with a CVODE solver
        boost::shared_ptr<AbstractCvodeCell> p_model(new Celldokos_model_1996FromCellMLCvode(p_solver, p_stimulus));

        std::cout << "\nModel = " << p_model->GetSystemName() << "\n";

        // Get the control funny current conductances for sodium and potassium, IK and ICaL.
        const double gf_Na = p_model->GetParameter("membrane_hyperpolarisation_activated_funny_current_sodium_component_conductance");
        const double gf_K = p_model->GetParameter("membrane_hyperpolarisation_activated_funny_current_potassium_component_conductance");
        const double g_K = p_model->GetParameter("membrane_delayed_rectifier_potassium_current_conductance");
        const double g_CaL = p_model->GetParameter("membrane_L_type_calcium_current_conductance");

        // Now we'll do four experiments for fun - block of each current, and then block of them in the correct proportions.
        for (unsigned channels=0u; channels<4u; channels++)
        {
            std::cout << "\nNew experiment: " << channels  << "\n";
            for (unsigned block_on=0; block_on<2u; block_on ++)
            {
                // Reset to defaults
                p_model->SetParameter("membrane_hyperpolarisation_activated_funny_current_sodium_component_conductance", gf_Na);
                p_model->SetParameter("membrane_hyperpolarisation_activated_funny_current_potassium_component_conductance", gf_K);
                p_model->SetParameter("membrane_delayed_rectifier_potassium_current_conductance", g_K);
                p_model->SetParameter("membrane_L_type_calcium_current_conductance", g_CaL);

                if (block_on)
                {
                    if (channels==0u || channels==3u)
                    {
                        double scaling = 0.81;
                        p_model->SetParameter("membrane_hyperpolarisation_activated_funny_current_sodium_component_conductance", gf_Na*scaling);
                        p_model->SetParameter("membrane_hyperpolarisation_activated_funny_current_potassium_component_conductance", gf_K*scaling);
                    }

                    if (channels==1 || channels==3u)
                    {
                        double scaling = 0.65;
                        p_model->SetParameter("membrane_delayed_rectifier_potassium_current_conductance", g_K*scaling);
                    }

                    if (channels==2 || channels==3u)
                    {
                        double scaling = 0.88;
                        p_model->SetParameter("membrane_L_type_calcium_current_conductance", g_CaL*scaling);
                    }
                }

                c_vector<double,2> results = RunToSteadyStateGetCycleLengthAndApd("3uM", p_model, block_on, channels);
                std::cout << "Block on = " << block_on << ",\tCycle length = " << results[1] << "ms,\tAPD50 = " << results[0] << "ms\n";
            }

        }
    }

    /**
     * Calculate a scaling for the maximum conductance, based on
     * IC50 for a channel and a requested concentration of a compound.
     *
     * (assumes Hill slope = 1)
     *
     * @param ic50  The IC50 value for this channel
     * @param conc  The concentration we are handling
     */
    double CalculateScaling(double ic50, double conc)
    {
        return 1.0 / (1.0 + (conc/ic50)/* to the power of Hill slope==1 in this case*/);
    }

};

#endif // TESTSANWITHFUNNYCURRENTBLOCK_HPP_
