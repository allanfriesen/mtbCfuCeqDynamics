#This code defines a first order ODE kinetic model to be used as a basis for SSA simulations
class KineticModel:

    def __init__(self, name='model'):
        self.name = name
        self.initialPopulations = []
        self.v = [] #v = state change vector
        self.rates=[]
        self.reactants=[]
        self.switchTimes = []
        self.speciesNames = []
        self.nSpecies = 0
        self.nRxn = 0
        self.nStages = 0

    def add_species(self, speciesName):
        self.nSpecies += 1
        self.speciesNames.append(speciesName)

    def add_stage(self, stageRates, switchTime): #stage ENDS at switch time
        if (len(stageRates) != self.nRxn):
            exit(' '.join(('Rate vector has',str(len(stageRates)),'elements, but there are',str(self.nRxn),
                'reactions in the model. Do not add rate vectors until all species are added to the model.')))

        for s in range(self.nRxn):
            self.rates[s].append(stageRates[s])
        self.nStages += 1
        self.switchTimes.append(switchTime)

    def add_reaction(self,reactant,v):
        if (len(v) != self.nSpecies):
            exit(' '.join(('State change vector has',str(len(v)),'elements, but there are',str(self.nSpecies),
                'species in the model. Do not add reaction vectors until all species are added to the model.')))
        self.v.append(v)
        self.reactants.append(reactant)
        self.nRxn += 1
        self.rates.append([])

    def get_nSpecies(self):
        return self.nSpecies

    def get_speciesNames(self):
        return self.speciesNames

    def printSpecies(self):
        print("nSpecies = ", self.nSpecies)
        print(self.speciesNames)
