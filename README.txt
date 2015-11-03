%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CC Attribution-NonCommercial-ShareAlike 4.0 International 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function fits HH equations on mechanically activated currents in
% various neurons.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   Requirements

   Current implementation requires Matlab 2014b or newer. Limiting factor
   is parallel computing toolbox, because starting parallel pool requires
   different syntax in older Matlab versions (2011b - 2014b).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

   Quick instructions

   1)  Copy the function in desired folder.
   2)  One starts the function with entering "Modeling_MS('Test')" at the 
       command window workspace MUST be set to folder containing 
       Modeling_MS.m
   3)  Folder containing script must contain experimental subfolder (in 
       this example case 'Test') in which the two files reside (see 
       example):
       -file with start-up parameters (Test_InitialFitParameters.mat)
       -file with the initial fit parameters (Test_ExperimentalData.mat)
   4)  The experimental subfolder must also contain folder 'Results' in 
       which the program dumps the results, figures etc....

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Several tests were implemented in this model, following guide lines in
   excellent paper of Hao & Delmas, 2010, which contains probably biggest 
   set of experiments regarding desensitisation mechanisms.
   
   Model is designed to fit several types of experiments, singly or 
   simultaneously:
       1)  simple trapezoidal stimuli (ramp-and-hold-and release)
               1.1)    current traces evoked by the stimulus
               1.2)    I-R curve of the above traces
       2)  channel inactivation with a double pulse protocol. First pulse 
           of varying amplitude is used for conditioning. Test pulse is  
           delivered immediately after conditioning with a saturating 
           amplitude. Test reveals a fraction of channels inactivated 
           by conditioning.
               2.1)    inactivated fraction (peak currents amplitudes to
                       the test stimulus, plotted agains the conditioning
                       amplitude.
       3)  recovery from adaptation with a double pulse protocol. First
           pulse is delivered with saturating amplitude. Second, test
           stimulus is delivered with increasing time delay but with the 
           same amplitude.
               3.1)    current traces evoked by the stimulus
               3.2)    current peaks elicited by test stimulus plotted
                       against time delay
       4)  desensitisation with a double step protocol. Conditioning is
           done with a trapezoidal stimulus, on top of which a test
           stimulus of variable amplitude is added at various time
           points.
               4.1)    time constants of inactivation and adaptation
               4.2)    inactivated fraction
               4.3)    adaptive shift
   
   An investigator can use this model to fit her or his own data (or
   perform meta-analysis), but there are minimum data requirements (see
   below):

   Minimum model data requirements:
       experimental data:  current traces and stimulus description
                           intensity-response(I-R) curve from exp. data

   Switches, turning on/off fitting of the selected experiment are 
   described below, in the section concerning fit parameters (C_Weights).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Included data set was obtained by running our model with various
   parameter set.

   To test your own data, use the structure of test InitialFitParam and
   ExperimentalData.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   Fit parameters are stored in XX_InitialFitParameters.m 
       XX is replaced by project folder name
   
   1)    C_Weights:  a "switch" vector (length = 14) turning on or off fitting
           of experiments and rules. To turn off fitting replace 1 with 0.
           See description of experiments above. If your data does not
           contain certain experiments, you can still use the model by
           turning off fitting of selected experiments.

           Fitting simple trapezoidal stimuli (ramp-and-hold-and release) 
           and its I-R curve. Serves as a control.
               C_Weights(1) .. compare current traces          
               C_Weights(7) .. restrict current traces to 
                               approach 0 mV in steady-state
               C_Weights(9) .. keep current traces before the stimulus in 
                               above 0.02
               C_Weights(2) .. compare I-R curves
               C_Weights(10).. penalise if min peak goes below -1.1

           Fitting recovery from adaptation with a double pulse protocol
               C_Weights(3) .. compare current traces          
               C_Weights(12).. compare recovery peaks of intensity    
                               response curves

           Fitting channel inactivation with a double pulse protocol
               C_Weights(14).. compare the I-R curves between control 
                               (simple trapezoidal stimulus) and inactivated
                               with prepulse

           Fitting desensitisation with a double step protocol.
               C_Weights(4) .. measures difference between modelled 
                               Tau Inact & Tau Adapt and measured 
                               Tau Inact & Tau Adapt
               C_Weights(5) .. measures difference between modelled 
                               adaptive shift and measured adaptive shift
               C_Weights(6) .. measures difference between modelled 
                               inactivated fraction and measured inactivated
                               fraction

           Restricting activation and inactivation
               C_Weights(8) .. penalise if h0 goes above 0.1 at x = 0
               C_Weights(13).. penalise if h0 goes above 0.1 at x = 0

   2)  dt                   .. integration time step (determining sample 
                               rate)
   3)  ExperimentalData     .. name of the file containing data used in
                               fitting
   4)  Model                .. name of the model used in fitting (only one
                               implemented in the current version)
   5)  ModelInput_FitVariables ..vector containing initial values for
                               variables used in fitting
   6)  ModelInput_FitVariables_names ..vector containing names of the
                               variables used in fitting
   7)  ModelInput_ConstVariables ..vector containing values of the
                               constants used in fitting
   8)  ModelInput_ConstVariables_names ..vector containing names of the
                               constants used in fitting
   9)  ModelInput_FitVariables_limits ..vector containing Hi and Lo
                               limits which limit possible values a 
                               variable can assume
   10) ModelInput_FitVariables_tol ..vector containing tolerance values 
                               (in fractions) if variable goes over the 
                               limit
   11) variableFitC_Weights .. Executing in multiple loop: matrix
                               containing multiple switch vectors for
                               multiple iterations of fitting. N of
                               columns should match the number of desired
                               iterations.
   11) variableFitParams    .. Executing in multiple loop: matrix
                               containing multiple vectors of initial fit 
                               parameters for multiple iterations of 
                               fitting. N of columns should match the 
                               number of desired iterations.
   12) variableFitParams_HiLims.. Executing in multiple loop: matrix
                               containing multiple vectors of Hi limits 
                               for multiple iterations of fitting. N of
                               columns should match the number of desired
                               iterations.
   13) variableFitParams_LoLims.. Executing in multiple loop: matrix
                               containing multiple vectors of Lo limits 
                               for multiple iterations of fitting. N of
                               columns should match the number of desired
                               iterations.
   14) variableFitParams_to .. Executing in multiple loop: matrix
                               containing multiple vectors of tolerances 
                               for multiple iterations of fitting. N of
                               columns should match the number of desired
                               iterations.
   15) C1 ..                   structure containing subsets for 
                               fitting simple trapezoidal stimuli 
                               (switch C_Weight(1))
           StimulusTimeSpan .. determining the end of the relevance of the
                               data
           RecordingWeights .. weight values for individual stimulus value
           IndexOfStimuliToFit..switch turning on/off fitting of certain
                               stimuli
   16) C2 ..                   structure containing subsetting for fitting
                               I-R curve (switch C_Weight(2)
           PointWeights ..     weight values for individual peak currents
   17) C3 ..                   structure containing subsets for 
                               fitting current responses 
                               (switch C_Weight(3))
           StimulusTimeSpan .. determining the end of the relevance of the
                               data
           RecordingWeights .. weight values for individual stimulus value
           IndexOfStimuliToFit..switch turning on/off fitting of certain
                               stimuli
           PointWeights ..     weight values for individual peak currents
   18) C4 ..                   structure containing subsets for 
                               fitting time constants (switch C_Weight(4))
           IndexOfStimuliToFit..switch turning on/off fitting of certain
                               stimuli
           PointWeights ..     weight values for individual peak currents
   19) C14 ..                  structure containing subsets for
                               fitting figure "Channel inactivation with a 
                               double pulse protocol (switch C_Weight(14))
           StimulusTimeSpan .. determining the end of the relevance of the
                               data
           RecordingWeights .. weight values for individual stimulus value
           IndexOfStimuliToFit..switch turning on/off fitting of certain
                               stimuli

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Current implementation fits 8 variables:
       tau_m   ..  time constant of activation
       tau_h   ..  time constant of inactivation
       A       ..  scaling of the stimulus
       B       ..  factor describing the linear relationship between 
                   adaptation and inactivation
       Hm      ..  slope for activation Boltzmann
       Xm      ..  midpoint for activation Boltzmann
       Hh      ..  slope for inactivation Boltzmann
       Xh      ..  midpoint for inactivation Boltzmann

   Besides, 6 constants are used (only 4 are currently active)
       Am      ..  scaling for activation Boltzmann
       Ah      ..  scaling for inactivation Boltzmann
       N       ..  power reflecting the number of channel subunits
                   (inactivation)
       M       ..  power reflecting the number of channel subunits
                   (activation)

   Constants can be easily transferred to variables and vice versa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Data that will be fitted are stored in XX_ExpementalData.mat
   !!Name of the file can be changed in InitialFitParam.mat!! 
   Currently the file contains data set was obtained by running our model
   with various parameter sets.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
