from typing import List, Callable, Union, Any
from .ga import GaData, MultiGaData, BinaryPopulation, BinaryGrayPopulation, VEGAHyperbolaFitness, BasePenalty, \
    VEGATournamentSelection, VEGABalancedSelection, VEGARankedSelection, OnePointCrossover, BinStringMutation, \
    NewGeneration, MultiNewGeneration


class VEGA:
    """
    Class for perform VEGA algorithm for multiobjective optimization.
    """
    def __init__(self, num_generations: int, num_individ: int, gens: tuple, obj_function: Callable,
                 obj_value: Union[int, float] = None, input_data: Any = None, penalty: BasePenalty = None,
                 conditions: list = None, is_gray: bool = False, children_percent: float = 0.95,
                 early_stop: Union[int, None] = 10, input_population: List[list] = None, tournament_size: int = 5,
                 mutation_lvl: Union[str, float] = 'normal', transfer_parents: str = 'best',
                 is_print: bool = True) -> None:
        """
        Initialization parameters:
        :param num_generations: int, number of generations;
        :param num_individ: int, number of individuals in generation (size of population);
        :param gens: tuple, tuple of 3 integer value, example: (0, 7, 1), 0 - min value,
                     7 - maximum value, 1 - step;
        :param obj_function: Callable, object function with 1 or 2 arguments, my_func(gens: list) or
                             my_func(input_data: Any, gens: list);
        :param obj_value: int | float, default: None. If object value exists, GA will optimize to the value,
                          else GA will optimize to min;
        :param input_data: Any, default: None. Argument for object function;
        :param penalty: BasePenalty, default: None. Subclass of BasePenalty(), initialization before
                        initialization subclass of BaseFitness(), used for conditional optimization.
                        Example: DynamicPenalty();
        :param conditions: list of strings (optimizer and conditionals) 3 value can be use: 'optimize', '<=', '!='.
                           Default: None.
                           Example: There is objective function: my_obj_func(x1, x2):
                                                                    return x1**2 + x2**2, 1-x1+x2, x1+x2
                                    my_obj_func returns 3 values, first value to optimize, second value must be <= 0,
                                    third value != 0, so have conditions = ['optimize', '<=', '!=']
        :param is_gray: bool, default: False. Ability to use gray code instead of binary representation;
        :param children_percent: float, default: 0.95. Percent of children who will be in new generation;
        :param early_stop: int, default: 10. Early stopping criteria, number of generation without improve;
        :param input_population: list[list], default: None. First generation from user;
        :param tournament_size: int, default: 3. Size of tournament, use only with selection_type="tournament";
        :param mutation_lvl: str | float, default: 'normal'. Mutations gens with different parameters:
                             float value: 0.00,...,1.00;
                             str value: 'weak', 'normal', 'strong';
        :param transfer_parents: str, default: "best". Type of transfer parents: "best", "random";
        :param is_print: bool, default: True. If True printed best solution;
        :return None
        """
        self.num_generations = num_generations
        self.num_individ = num_individ
        self.gens = gens
        self.obj_function = obj_function
        self.obj_value = obj_value
        self.input_data = input_data
        self.penalty = penalty
        self.conditions = conditions
        self.is_gray = is_gray
        self.children_percent = children_percent
        self.early_stop = early_stop if early_stop is not None else num_generations
        self.input_population = input_population
        self.tournament_size = tournament_size
        self.mutation_lvl = mutation_lvl
        self.transfer_parents = transfer_parents
        self.is_print = is_print

    def optimize(self) -> GaData:
        """
        Main method of VEGA.
        :return: GaData
        """
        # init GaData & Population
        num_objectives = len(self.obj_function([0]*len(self.gens))) if self.conditions is None\
            else self.conditions.count('optimize')
        ga_data = MultiGaData(num_generations=self.num_generations, children_percent=self.children_percent,
                              early_stop=self.early_stop)
        if self.is_gray:
            population = BinaryGrayPopulation()
        else:
            population = BinaryPopulation()

        population.set_params(num_individ=self.num_individ, gens=self.gens, input_population=self.input_population)
        # init fitness func, selection, crossover, mutation, new generation
        fitness_func = VEGAHyperbolaFitness(obj_function=self.obj_function, obj_value=self.obj_value,
                                            input_data=self.input_data, penalty=self.penalty,
                                            conditions=self.conditions)
        selection = VEGATournamentSelection(num_objectives=num_objectives)
        cross = OnePointCrossover()
        mutation = BinStringMutation(mutation_lvl=self.mutation_lvl)
        new_generation = MultiNewGeneration(transfer_parents=self.transfer_parents)
        # creating first generation
        population.fill()
        ga_data.population = population
        fitness_func.execute(ga_data=ga_data)
        ga_data.update()
        # main loop for GA perform
        for i in range(1, ga_data.num_generations):
            selection.execute(ga_data)
            cross.execute(ga_data)
            mutation.execute(ga_data)
            new_generation.execute(ga_data)
            fitness_func.execute(ga_data)
            ga_data.update()

            if ga_data.num_generation_no_improve > ga_data.early_stop:
                break

        if self.is_print:
            ga_data.print_best_solution()
        return ga_data
