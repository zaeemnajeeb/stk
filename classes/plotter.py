import os
import matplotlib.pyplot as plt

from ..fitness import *

class Plotter:
    """
    A class for handling any and all plotting of a population's data.    
    
    Each population is initialized with a ``Plotter`` instances in its
    `plot` attribute. In turn, the plotter instance holds the population
    in its `pop` attribute.    
    
    Attributes
    ----------
    pop : Population
        The population whose data a ``Plotter`` instance plots.
    
    """

    def __init__(self, pop):
        self.pop = pop
    
    def epp(self, plot_name):

        progress = self.pop.ga_tools.ga_progress        
        
        if (isinstance(progress.maxs[0], float) or 
                              isinstance(progress.maxs[0], int)):
            fig = plt.figure()
            plt.xlabel('Generation Number')
            plt.ylabel('Fitness Value')
            plt.title('Evolutionary Progress Plot', fontsize=18)
            plt.plot(progress.gens, progress.means, color='green')
            plt.plot(progress.gens, progress.mins, color='blue')
            plt.plot(progress.gens, progress.maxs, color='red')
            fig.savefig(plot_name, dpi=1000)
            plt.close('all')
            
        else:
            func_name = self.pop.ga_tools.ga_input.fitness_func.name
            fitness_func = globals()[func_name]
            
            for x in range(len(progress.means[0])):
                fig = plt.figure()
                plt.xlabel('Generation Number')                             
                plt.ylabel('Unscaled ' + fitness_func.param_labels[x])                
                plt.title(' Evolutionary Progress Plot', fontsize=18)
                y_mean = [v[x] for v in progress.means]
                y_max = [v[x] for v in progress.maxs]
                y_min = [v[x] for v in progress.mins]

                plt.plot(progress.gens, y_mean, color='green')
                plt.plot(progress.gens, y_min, color='blue')
                plt.plot(progress.gens, y_max, color='red')
                
                new_plot_name = str(x).join(os.path.splitext(plot_name))                
                
                fig.savefig(new_plot_name, dpi=1000)
                plt.close('all')

    def subpopulations(self, plot_name):
        """
        Plots the min, max and avg fitness values of each subpopulation.
        
        Parameters
        ----------
        plot_name : str
            The full path of where the plot should be saved.
            
        Returns
        -------
        None : NoneType
        
        """
        
        xvals = []
        maxs = []
        means = []
        mins = []
        
        for i, subpop in enumerate(self.pop.populations, 1):
            xvals.append(i)
            maxs.append(max(x.fitness for x in subpop))
            means.append(subpop.mean(lambda x : x.fitness))
            mins.append(min(x.fitness for x in subpop))
        
        fig = plt.figure()
        plt.xlabel('Population')
        plt.ylabel('Fitness')
        plt.scatter(xvals, maxs, color='red', marker='x')
        plt.scatter(xvals, means, color='green', marker='x')
        plt.scatter(xvals, mins, color='blue', marker='x')
        fig.savefig(plot_name, dpi=1000)
        plt.close('all')
    
    def progress_params(self, plot_name):
        """
        Plots the progress_params values across subpopulations.

        For each progress_param a plot will be produced. Each will show
        subpopulations on the x-axis and the min, max and avg values
        of that progress param on the y axis.

        Parameters
        ----------
        plot_name : str
            The full path of where the plots should be saved.
            
        Returns
        -------
        None : NoneType
        
        """
        
        progress = self.pop.ga_tools.ga_progress        
        
        func_name = self.pop.ga_tools.ga_input.fitness_func.name
        fitness_func = globals()[func_name]
        
        min_params = []
        max_params = []
        mean_params = []
        
        for subpop in self.pop.populations:
            min_params.append(subpop.progress.mins[-1])
            max_params.append(subpop.progress.maxs[-1])
            mean_params.append(subpop.progress.means[-1])
        
        for x in range(len(params[0])):
            fig = plt.figure()
            plt.xlabel('Generation Number')                             
            plt.ylabel('Unscaled ' + fitness_func.param_labels[x])                
            plt.title('Population Comparison', fontsize=18)
            y_mean = [array[x] for array in mean_params]
            y_max = [array[x] for array in max_params]
            y_min = [array[x] for array in min_params]

            plt.scatter(progress.gens, y_mean, color='green', marker='x')
            plt.scatter(progress.gens, y_min, color='blue', marker='x')
            plt.scatter(progress.gens, y_max, color='red', marker='x')
            
            new_plot_name = str(x).join(os.path.splitext(plot_name))                
            
            fig.savefig(new_plot_name, dpi=1000)
            plt.close('all')
        
        
        