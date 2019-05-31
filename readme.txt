This is a quick start guide for how to use iGrow, a compartmental
simulation of diffusion and active transport in a neuron. You need
python with numpy installed.

Run the examples from our article:


Figure 1:

Illustration of model, no simulations to run.


Figure 2: 

Run the following matlab script to generate figure 2B1 and 2B2. The
script writes a list of jobs to the file
FIG2A-diffusion-only-summary.txt located in the "input" directory. You
can in the script specify how many workers you want to have perform
the jobs. The default is nWorkers = 3. If you increase this number
then you need to start more workers to make sure that all the jobs get
run as each job belongs to a specific worker. Running the setup script
creates a neuron with one primary branch and two secondary branches,
and repeats the simulations with different diffusion constants.

runSimFig2A.m

The script will display three lines of commands that are used to run
the python workers that do the simulations:

python runSimWorker.py input/FIG2A-diffusion-only-summary.txt 1
python runSimWorker.py input/FIG2A-diffusion-only-summary.txt 2
python runSimWorker.py input/FIG2A-diffusion-only-summary.txt 3

After the simulations are done run the following script to generate
the figures 2B1 and 2B2.

plotTwoGrowthConesFig2A.m

Status: Verified.


To generate figures 2C1 and 2C2 please run the following script, and
then start the python workers indicated (as per the example above):

runSimFig2C.m

To run workers:

python runSimWorker.py input/Fig2C-diffusion-and-active-transport-summary.txt 1
python runSimWorker.py input/Fig2C-diffusion-and-active-transport-summary.txt 2
python runSimWorker.py input/Fig2C-diffusion-and-active-transport-summary.txt 3

Generate the figures for Figure 2C1 and 2C2:

plotTwoGrowthConesFig2C.m


Finally to make the Figure 2D, start by setting up the work list:

runSimFig2D.m

Run the workers:

python runSimWorker.py input/Fig2D-X-Y-range-summary.txt 1
python runSimWorker.py input/Fig2D-X-Y-range-summary.txt 2
python runSimWorker.py input/Fig2D-X-Y-range-summary.txt 3

Generate the figures:

plotTwoGrowthCones2D.m


Status: Verified.

Figure 3:

This figure shows results of competition in a model simulation
starting from a real hippocampal neuron morphology. The simulation is
run as many times as there are growth cones, plus one reference
case. In the line below, replace the data path with your own path to
the files, also the number at the end must be changed to run different
growth cones (-1 is reference).

python runSimDistanceDependence.py DATA/cell1.swc 0

To generate the dendrogram figure, copy the data files, and run:

analyseDistanceDependenceCompetition.m

To get the rendering of the morphology you need to have povray installed.

exportToPovRay('DATA/','cell1.swc-NOGCmod-out.txt',[],500)

Status: Partly verified, no povray on my mac machine.


Figure 4:

To initialise the simulations, run the following matlab scripts:

setupPredictSpeedSimulation.m

You have the option to use either the standard workers, or try the
parallell version of the script:

python runSimWorkerPredictSpeed.py input/GridSearch/predict-growth-speed-summary.txt 1

python runSimWorkerPredictSpeedParallel.py input/GridSearch/predict-growth-speed-summary.txt

The result figures are generated using:

calculatePredictionFitness

Status: Partly verified, ran a limited search.


