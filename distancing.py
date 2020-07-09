#! /usr/bin/env python3

from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path
import json
import matplotlib.pyplot as plt
import numpy
import os
import gillespy2
import stochpy
import tempfile

from gillespy2.solvers import NumPySSASolver
from gillespy2.solvers import TauHybridSolver


R0 = 2.4
MEAN_INCUBATION_PERIOD = 4
SYMPTOMATIC_FRACTION = 1/5
MEAN_INFECTIOUS_PERIOD = 6

TOTAL_POPULATION = 10000
INITIAL_INFECTIOUS = 2
POPULATIONS = {
    "S": TOTAL_POPULATION-INITIAL_INFECTIOUS,
    "E": 0,
    "Y": 0,
    "A": 0,
    "I": INITIAL_INFECTIOUS,
    "C": INITIAL_INFECTIOUS,
    "R": 0
}

START_TIME = 0
END_TIME = 180
TIMESTEPS = (END_TIME - START_TIME) + 1
REALIZATIONS = 10


def main(do_cms, do_gillespy, do_stochpy):

    if do_cms:
        run_cms()

    if do_gillespy:
        run_gillespy()

    if do_stochpy:
        run_stochpy()

    return


def run_cms():

    emodl = f"""
    (import (rnrs) (emodl cmslib))
    (start-model "seir")
    (species S {POPULATIONS["S"]})
    (species E {POPULATIONS["E"]})
    (species Y {POPULATIONS["Y"]})
    (species A {POPULATIONS["A"]})
    (species I {POPULATIONS["I"]})
    (species C {POPULATIONS["C"]})
    (species R {POPULATIONS["R"]})
    (param Ki {R0/MEAN_INFECTIOUS_PERIOD})
    (param Ke {1/MEAN_INCUBATION_PERIOD})
    (param Ky {SYMPTOMATIC_FRACTION/MEAN_INCUBATION_PERIOD})
    (param Ka {(1-SYMPTOMATIC_FRACTION)/MEAN_INCUBATION_PERIOD})
    (param Ks {SYMPTOMATIC_FRACTION})
    (param Kr {1/MEAN_INFECTIOUS_PERIOD})
    (param V 0)
    (reaction transmission (S) (E)     (/ (* Ki S I) (+ S E I R)))
    (reaction infectious_a (E) (I A C) (* Ka E))
    (reaction infectious_s (E) (I Y C) (* Ky E))
    (reaction recovery     (I) (R)     (* Kr I))
    (state-event sia (> Y 100) ((Ki {0.8/MEAN_INFECTIOUS_PERIOD})))
    (observe susceptible  S)
    (observe exposed      E)
    (observe symptomatic  Y)
    (observe asymptomatic A)
    (observe infectious   I)
    (observe cases        C)
    (observe recovered    R)
    (end-model)
    """
    t0 = datetime.now()
    model = EmodlLoader.LoadEMODLModel(str(emodl))

    config = {
        "solver": "SSA",
        "runs": REALIZATIONS,
        "duration": END_TIME-START_TIME,
        "samples": TIMESTEPS,
        "prng_seed": datetime.now().microsecond
    }
    cfg.CurrentConfiguration = cfg.ConfigurationFromString(json.dumps(config))

    solver = solvers.CreateSolver(config["solver"], model, config["runs"], config["duration"], config["samples"])
    t1 = datetime.now()
    """
    Solve() is a wrapper for:
        for (int curRealization = 0; curRealization < numRealizations; curRealization++)
        {
            StartRealization(); // sets CurrentTime to 0.0 and calls ResetModelState()
            SolveOnce();
        }
        
    SolveOnce() is a wrapper for:
        while (CurrentTime < duration)
        {
            StepOnce();
        }
        
    StepOnce() looks like this:

        double timeOfNextEvent = ExecuteScheduledEvents(duration);
        double newTau = CalculateProposedTau(timeOfNextEvent);
        CurrentTime = newTau;
        SamplingParams = trajectories.RecordObservables(model.Observables, SamplingParams, CurrentTime, duration);
        if (CurrentTime < duration) 
        {
            ExecuteReactions();
            UpdateTriggeredEvents();
        }
    """

    solver.Solve()
    t2 = datetime.now()

    print(f"Time for model construction: {t1-t0}")
    print(f"Time for model execution ({REALIZATIONS} trajectories): {t2-t1}")

    data = solver.GetTrajectoryData()
    for index, label in enumerate(solver.GetTrajectoryLabels()):
        plt.plot([float(value) for value in data[index]], label=str(label))

    plt.title("PyCMS")
    plt.legend()
    plt.show()

    return


class SEIR(gillespy2.Model):

    def __init__(self):
        super().__init__(self)
        s = gillespy2.Species(name="susceptible", initial_value=POPULATIONS["S"])
        e = gillespy2.Species(name="exposed", initial_value=POPULATIONS["E"])
        y = gillespy2.Species(name="symptomatic", initial_value=POPULATIONS["Y"])
        a = gillespy2.Species(name="asymptomatic", initial_value=POPULATIONS["A"])
        i = gillespy2.Species(name="infectious", initial_value=POPULATIONS["I"])
        c = gillespy2.Species(name="cases", initial_value=POPULATIONS["C"])
        r = gillespy2.Species(name="recovered", initial_value=POPULATIONS["R"])
        self.add_species([s, e, y, a, i, c, r])
        k_i = gillespy2.Parameter(name="Ki", expression=R0/MEAN_INFECTIOUS_PERIOD)
        # k_e = gillespy2.Parameter(name="Ke", expression=1/MEAN_INCUBATION_PERIOD)
        k_y = gillespy2.Parameter(name="Ky", expression=SYMPTOMATIC_FRACTION/MEAN_INCUBATION_PERIOD)
        k_a = gillespy2.Parameter(name="Ka", expression=(1-SYMPTOMATIC_FRACTION)/MEAN_INCUBATION_PERIOD)
        k_s = gillespy2.Parameter(name="Ks", expression=SYMPTOMATIC_FRACTION)
        k_r = gillespy2.Parameter(name="Kr", expression=1/MEAN_INFECTIOUS_PERIOD)
        vaccinated = gillespy2.Parameter(name="vaccinated", expression=0)
        # self.add_parameter([k_i, k_e, k_s, k_r])
        self.add_parameter([k_i, k_y, k_a, k_s, k_r, vaccinated])
        transmission = gillespy2.Reaction(name="transmission",
                                          propensity_function="Ki*susceptible*infectious/(susceptible+exposed+infectious+recovered)",
                                          reactants={s: 1}, products={e: 1})
        transmissive_a = gillespy2.Reaction(name="transmissive_asymptomatic",
                                            massaction=True, rate=k_a,
                                            reactants={e: 1}, products={a: 1, i: 1, c: 1})
        transmissive_s = gillespy2.Reaction(name="transmissive_symptomatic",
                                            massaction=True, rate=k_y,
                                            reactants={e: 1}, products={y: 1, i: 1, c: 1})
        recovery = gillespy2.Reaction(name="recovery", massaction=True, rate=k_r,
                                      reactants={i: 1}, products={r: 1})
        self.add_reaction([transmission, transmissive_a, transmissive_s, recovery])

        initiate_distancing = gillespy2.EventAssignment(variable=k_i, expression=f"{0.8/MEAN_INFECTIOUS_PERIOD}")
        trigger = gillespy2.EventTrigger(expression="symptomatic > 100", initial_value=False, persistent=True)
        sia = gillespy2.Event("distancing", delay="0", assignments=[initiate_distancing],
                              trigger=trigger, use_values_from_trigger_time=False)
        self.add_event([sia])

        self.timespan(numpy.linspace(0, 180, 181))

        return


def run_gillespy():

    t0 = datetime.now()
    model = SEIR()
    t1 = datetime.now()
    # NumPySSASolver does not support events
    results = model.run(solver=TauHybridSolver, number_of_trajectories=REALIZATIONS)
    t2 = datetime.now()

    print(f"Time for model construction: {t1-t0}")
    print(f"Time for model execution ({REALIZATIONS} trajectories): {t2-t1}")

    for index in range(REALIZATIONS):
        trajectory = results[index]
        plt.plot(trajectory["time"], trajectory["susceptible"], label="susceptible")
        plt.plot(trajectory["time"], trajectory["exposed"], label="exposed")
        plt.plot(trajectory["time"], trajectory["infectious"], label="infectious")
        plt.plot(trajectory["time"], trajectory["recovered"], label="recovered")
        plt.plot(trajectory["time"], trajectory["symptomatic"], label="symptomatic")
        plt.plot(trajectory["time"], trajectory["asymptomatic"], label="asymptomatic")
        plt.plot(trajectory["time"], trajectory["cases"], label="total cases")

    plt.title("GillesPy2")
    plt.legend()
    plt.show()

    return


def run_stochpy():

    pysces = """
    Modelname: stochpy_seir
    Description: SEIR model with COVID-19 like parameters

    # Symbol names (i.e. reaction, species, compartment, function, rule and parameter names etc.) must start with either
    # an underscore or letter and be followed by any combination of alphanumeric characters or an underscore. Like all
    # other elements of the input file names are case sensitive.
    # Explicit access to the "current" time in a time simulation is provided by the special symbol _TIME_.

    # Compartment: <name>, <size>,
    # <dimensions>, where <name> is the unique compartment id, <size> is the size of the
    # compartment (i.e. length, volume or area) defined by the number of <dimensions> (e.g. 1,2,3)

    # Function: <name>, <args> {
    #     <formula>
    # }
    # where <name> is the unique function id, <arglist> is one or more comma separated function arguments. The <formula>
    # field, enclosed in curly brackets, may only make use of arguments listed in the <arglist> and therefore cannot
    # reference model attributes directly. If this functionality is required a forcing function (assignment rule) may be
    # what you are looking for.

    # The reaction stoichiometry and rate equation are defined together as a single reaction step. Each step in the
    # system is defined as having a name (identifier), a stoichiometry (substrates are converted to products) and rate
    # equation (the catalytic activity, described in terms of species and parameters).
    # The PySCeS MDL also allows the use of the $pool token that represents a placeholder reagent for reactions that
    # have no net substrate or product.
    # <name>[@<compartment>]:
    #     <stoichiometry> e.g. {2}S + I > S + 2{I}
    #     <rate equation>

    transmission:
        S > E
        Ki * S * I / (S + E + I + R)
        
    transmissive_a:
        E > I + A + C
        Ka * E
        
    transmissive_s:
        E > I + Y + C
        Ky * E
        
    recovery:
        I > R
        Kr * I

    # Standard Python operators + - * / ** are supported
    # the following functions are supported in any mathematical expression:
    # log, log10, ln, abs, pow, exp, root, sqrt, sin, cos, tan, sinh, cosh, tanh, floor, ceil, ceiling, piecewise, ...
    # Logical operators are supported in rules, events etc but _not_ in rate equation definitions.
    
    # The general form of any species (fixed, free) and parameter is simply: property = value
    
    S = 9998
    E =    0
    Y =    0
    A =    0
    I =    2
    C =    2
    R =    0
    
    Ki = 0.4   # 2.4 / 6
    Ke = 0.25  # 1.0 / 4
    Ka = 0.2   # 0.8 * Ke
    Ky = 0.05  # 0.2 * Ke
    Kr = 0.167 # 1.0 / 6
    
    # Assignment rules can access other model attributes directly and have the generic form !F <par> = <formula>.
    
    # The general form of an event is
    # Event: <name>, <trigger>, <delay> {
    #     <assignments>
    # }
    # As can be seen an event consists of essentially three parts, a conditional <trigger>, a set of one or more
    # <assignments> and a <delay> between when the trigger is fired (and the assignments are evaluated) and the eventual
    # assignment to the model. Assignments have the general form <par> = <formula>. Events have access to the
    # "current" simulation time using the _TIME_ symbol

    Event: sia, Y > 100, 0 {
        Ki = 0.32
    }

"""

    t0 = datetime.now()
    handle, filename = tempfile.mkstemp(suffix=".psc", text=True)
    os.write(handle, pysces.encode())
    os.close(handle)
    path = Path(filename).resolve()

    try:
        model = stochpy.SSA(model_file=path.name, dir=str(path.parent), end=END_TIME, trajectories=REALIZATIONS)
        t1 = datetime.now()
        model.DoStochSim(end=180, mode='time', method='tauleap', trajectories=REALIZATIONS, epsilon=0.025)
        t2 = datetime.now()

        print(f"Time for model construction: {t1 - t0}")
        print(f"Time for model execution ({REALIZATIONS} trajectories): {t2 - t1}")

        for t in range(REALIZATIONS):
            model.GetTrajectoryData(t+1)
            data = model.data_stochsim
            for i in range(data.species.shape[1]):
                plt.plot(data.time, data.species[:, i], label=data.species_labels[i])

        plt.title("StochPy")
        plt.legend()
        plt.show()
    finally:
        path.unlink()

    return


if __name__ == "__main__":

    CMS_PATH = os.environ["CMS_PATH"] if "CMS_PATH" in os.environ else None

    parser = ArgumentParser()
    parser.add_argument("-c", "--cms", default=True, action="store_false")
    parser.add_argument("-g", "--gillespy", default=True, action="store_false")
    parser.add_argument("-s", "--stochpy", default=True, action="store_false")
    parser.add_argument("-b", "--binary", default=CMS_PATH, help="Path containing 'compartments.exe'")

    args = parser.parse_args()

    if args.cms:
        if args.binary is None:
            raise RuntimeError("Location of CMS binary must be specified to run CMS model.")
        else:
            import clr
            clr.AddReference(args.binary+"/compartments.exe")
            global EmodlLoader, cfg, solvers
            from compartments.emodl import EmodlLoader
            from compartments import Configuration as cfg
            from compartments.emod.utils import SolverFactory as solvers

    main(args.cms, args.gillespy, args.stochpy)
