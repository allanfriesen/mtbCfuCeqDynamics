import KineticModel as km
import subprocess

class StochasticSimulation:

    infile="/home/allanfriesen/local/stochastic/stochastic.in" #on work machine
    #infile="/Users/allanfriesen/local/stochastic/stochastic.in" #on laptop

    #cmd="/home/allanfriesen/bin/stochastic.exe"
    cmd="/Users/allanfriesen/bin/stochastic.exe"
    def __init__(self, model, nTraj=1, tSim=100, dtRecip = 1, initialPopulations=None, mode='normal'):
        self.model = model
        self.nTraj = nTraj
        self.tSim = tSim
        self.dtRecip = dtRecip
        self.initialPopulations = initialPopulations
        self.trajectories = []
        self.times = []
        self.mode = mode

    def write_infile(self):
        
        def datPrint( vec ):
            f.write(" ".join(map(str,vec)))
            f.write("\n")

        with open(self.infile, "w") as f:
            f.write(str(self.nTraj))
            f.write("\n")
            f.write(str(self.tSim))
            f.write("\n")
            f.write(str(self.dtRecip))
            f.write("\n")
            f.write(str(self.model.nSpecies))
            f.write("\n")
            if( self.initialPopulations == None ): #default initial populations
                popVec = [1]
                popVec = popVec + [0] * (self.model.nSpecies - 1)
            for i in range(self.nTraj):
                datPrint(self.initialPopulations[i])
            f.write(str(self.model.nRxn))
            f.write("\n")
            f.write(str(self.model.nStages))
            f.write("\n")
            for i in range(self.model.nRxn):
                for s in range(self.model.nSpecies):
                    f.write( str(self.model.v[i][s]))
                    f.write("\n")
                f.write(str(self.model.reactants[i]))
                f.write("\n")
            for i in range(self.model.nRxn):
                datPrint( self.model.rates[i] )
            datPrint( self.model.switchTimes )

    def run_simulation(self):
        self.trajectories = []
        self.write_infile()
        lines = subprocess.check_output(self.cmd).splitlines()
        if( self.mode == 'debug' ):
            print( lines )
        lineList = []
        for line in lines:
            lineList.append(  [float(i) for i in line.decode('(utf8)').split()]  )
        self.times = lineList[0]
        line = 1
        for i in range(self.model.nSpecies):
            self.trajectories.append([])
        #self.trajectories = [[]]*self.model.nSpecies
        #print()
        #print("lines")
        #print(lines)
        #print()
        for i in range(self.nTraj):
            for s in range(self.model.nSpecies):
                self.trajectories[s].append(lineList[line])
                #print('appended line ', line, ' for species ', s, ' for simulation ', i)
                #print(self.trajectories)
                line += 1
        #print()
        #print()
        #print( 'total lines: ', len(lineList) )
        return {'species': self.model.speciesNames, 'times': self.times, 'trajectories': self.trajectories}
