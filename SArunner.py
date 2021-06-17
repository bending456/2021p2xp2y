import numpy as np 
import scipy
import pickle as pk 
import math 
import ScriptRunner as SR 

def ATPtest(changingVar,ODEfilename,time):
    counter = 0
    for i in np.arange(len(changingVar)):
        data   = SR.gotranMicroglia(sim_time      = time,
                                    ATP           = changingVar[i],
                                    output_name   = 'test1',
                                    ode_file_name = ODEfilename,
                                    data_name2    = 'mRNA_TNF',
                                    data_name3    = 'pAkt',
                                    data_name4    = 'S2',
                                    data_name5    = 'TNFae',
                                    removePickle  = 1,
                                    timePrint     = 0)
                    
        if counter == 0:
            dura     = data[0]
            Ca       = data[1]
            mRNATNF  = data[2]
            pakt     = data[3]
            dist     = data[4]
            TNFae    = data[5]
        else:
            Ca       = np.vstack([Ca,data[1]])
            mRNATNF  = np.vstack([mRNATNF,data[2]])
            pakt     = np.vstack([pakt,data[3]])
            dist     = np.vstack([dist,data[4]])
            TNFae    = np.vstack([TNFae,data[5]])

        counter = counter + 1

    return dura, Ca, mRNATNF, pakt, dist, TNFae 


def P2X4test(atp,changingVar,ODEfilename,time):
    counter = 0
    for i in np.arange(len(changingVar)):
        data   = SR.gotranMicroglia(sim_time      = time,
                                    ATP           = atp, 
                                    rhop2x4       = changingVar[i],
                                    output_name   = 'test1',
                                    ode_file_name = ODEfilename,
                                    data_name2    = 'mRNA_TNF',
                                    data_name3    = 'pAkt',
                                    data_name4    = 'S2',
                                    data_name5    = 'TNFae',
                                    removePickle  = 1,
                                    timePrint     = 0)
                    
        if counter == 0:
            dura     = data[0]
            Ca       = data[1]
            mRNATNF  = data[2]
            pakt     = data[3]
            dist     = data[4]
            TNFae    = data[5]
        else:
            Ca       = np.vstack([Ca,data[1]])
            mRNATNF  = np.vstack([mRNATNF,data[2]])
            pakt     = np.vstack([pakt,data[3]])
            dist     = np.vstack([dist,data[4]])
            TNFae    = np.vstack([TNFae,data[5]])

        counter = counter + 1 

    return dura, Ca, mRNATNF, pakt, dist, TNFae 


def P2X7test(atp,changingVar,ODEfilename,time):
    counter = 0
    for i in np.arange(len(changingVar)):
        data   = SR.gotranMicroglia(sim_time      = time,
                                    ATP           = atp,
                                    rhop2x7       = changingVar[i],
                                    output_name   = 'test1',
                                    ode_file_name = ODEfilename,
                                    data_name2    = 'mRNA_TNF',
                                    data_name3    = 'pAkt',
                                    data_name4    = 'S2',
                                    data_name5    = 'TNFae',
                                    removePickle  = 1,
                                    timePrint     = 0)
                    
        if counter == 0:
            dura     = data[0]
            Ca       = data[1]
            mRNATNF  = data[2]
            pakt     = data[3]
            dist     = data[4]
            TNFae    = data[5]
        else:
            Ca       = np.vstack([Ca,data[1]])
            mRNATNF  = np.vstack([mRNATNF,data[2]])
            pakt     = np.vstack([pakt,data[3]])
            dist     = np.vstack([dist,data[4]])
            TNFae    = np.vstack([TNFae,data[5]])

        counter = counter + 1 

    return dura, Ca, mRNATNF, pakt, dist, TNFae 

def P2Xtest(atp,changingVar,ODEfilename,time):
    counter = 0
    for i in np.arange(len(changingVar)):
        data   = SR.gotranMicroglia(sim_time      = time,
                                    ATP           = atp,
                                    rhop2x4       = changingVar[i],
                                    rhop2x7       = changingVar[i],
                                    output_name   = 'test1',
                                    ode_file_name = ODEfilename,
                                    data_name2    = 'mRNA_TNF',
                                    data_name3    = 'pAkt',
                                    data_name4    = 'S2',
                                    data_name5    = 'TNFae',
                                    removePickle  = 1,
                                    timePrint     = 0)
                    
        if counter == 0:
            dura     = data[0]
            Ca       = data[1]
            mRNATNF  = data[2]
            pakt     = data[3]
            dist     = data[4]
            TNFae    = data[5]
        else:
            Ca       = np.vstack([Ca,data[1]])
            mRNATNF  = np.vstack([mRNATNF,data[2]])
            pakt     = np.vstack([pakt,data[3]])
            dist     = np.vstack([dist,data[4]])
            TNFae    = np.vstack([TNFae,data[5]])

        counter = counter + 1 

    return dura, Ca, mRNATNF, pakt, dist, TNFae 

def P2Ytest(atp,changingVar,ODEfilename,time):
    counter = 0
    for i in np.arange(len(changingVar)):
        data   = SR.gotranMicroglia(sim_time      = time,
                                    ATP           = atp,
                                    rhop2yc        = changingVar[i],
                                    output_name   = 'test1',
                                    ode_file_name = ODEfilename,
                                    data_name2    = 'mRNA_TNF',
                                    data_name3    = 'pAkt',
                                    data_name4    = 'S2',
                                    data_name5    = 'TNFae',
                                    removePickle  = 1,
                                    timePrint     = 0)
                    
        if counter == 0:
            dura     = data[0]
            Ca       = data[1]
            mRNATNF  = data[2]
            pakt     = data[3]
            dist     = data[4]
            TNFae    = data[5]
        else:
            Ca       = np.vstack([Ca,data[1]])
            mRNATNF  = np.vstack([mRNATNF,data[2]])
            pakt     = np.vstack([pakt,data[3]])
            dist     = np.vstack([dist,data[4]])
            TNFae    = np.vstack([TNFae,data[5]])

        counter = counter + 1 

    return dura, Ca, mRNATNF, pakt, dist, TNFae 
